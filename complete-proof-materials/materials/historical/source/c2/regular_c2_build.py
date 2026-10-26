#!/usr/bin/env python3
"""Deterministically regenerate the exact regular-C2 LOCAL extraction.

Python 3.10+ standard library only. No numerical tolerances or inherited receipts.
Run beside the generated JSON files; existing generated files are replaced only
inside this newly named reconstruction directory. This is not a global lift.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json

ROOT = Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit('Run without -O: exact assertions are required.')


class Poly:
    def __init__(self, n, terms=None):
        self.n = n
        self.d = {tuple(e): F(c) for e, c in (terms or {}).items() if c}

    @classmethod
    def const(cls, n, x):
        return cls(n, {(0,) * n: F(x)})

    @classmethod
    def var(cls, n, i):
        e = [0] * n
        e[i] = 1
        return cls(n, {tuple(e): F(1)})

    def coerce(self, q):
        return q if isinstance(q, Poly) else Poly.const(self.n, q)

    def __add__(self, q):
        q = self.coerce(q)
        assert self.n == q.n
        d = dict(self.d)
        for e, c in q.d.items():
            d[e] = d.get(e, F(0)) + c
            if not d[e]:
                del d[e]
        return Poly(self.n, d)

    __radd__ = __add__

    def __neg__(self):
        return Poly(self.n, {e: -c for e, c in self.d.items()})

    def __sub__(self, q):
        return self + (-self.coerce(q))

    def __rsub__(self, q):
        return self.coerce(q) - self

    def __mul__(self, q):
        if not isinstance(q, Poly):
            q = F(q)
            return Poly(self.n, {e: q*c for e, c in self.d.items()})
        assert self.n == q.n
        d = {}
        for e, c in self.d.items():
            for f, b in q.d.items():
                g = tuple(x+y for x, y in zip(e, f))
                d[g] = d.get(g, F(0)) + c*b
        return Poly(self.n, d)

    __rmul__ = __mul__

    def __truediv__(self, q):
        return self * (F(1)/q)

    def __pow__(self, k):
        assert isinstance(k, int) and k >= 0
        out = Poly.const(self.n, 1)
        for _ in range(k):
            out = out * self
        return out

    def derivative(self, i):
        d = {}
        for e, c in self.d.items():
            if e[i]:
                f = list(e)
                f[i] -= 1
                d[tuple(f)] = c*e[i]
        return Poly(self.n, d)

    def substitute(self, args):
        assert len(args) == self.n
        m = args[0].n
        maxes = [max((e[i] for e in self.d), default=0) for i in range(self.n)]
        powers = [[a**k for k in range(maxes[i]+1)] for i, a in enumerate(args)]
        out = Poly(m)
        for e, c in self.d.items():
            term = Poly.const(m, c)
            for i, k in enumerate(e):
                if k:
                    term = term * powers[i][k]
            out = out + term
        return out

    def constant(self):
        return self.d.get((0,)*self.n, F(0))

    def l1_nonconstant(self):
        return sum((abs(c) for e, c in self.d.items() if any(e)), F(0))

    def records(self):
        return [[list(e), str(c)] for e, c in sorted(self.d.items())]

    def __eq__(self, other):
        return isinstance(other, Poly) and self.n == other.n and self.d == other.d


SOURCE_ORDER = ['d', 'n', 'm', 'C', 'A', 'D', 'alpha']
CONTACT_ORDER = ['u', 'v', 'W', 'X', 'Y']
CSTAR = tuple(map(F, ['7/8', '9/8', '1/18', '11/8', '9/8', '1/18']))
THETASTAR = tuple(map(F, ['1', '1', '1', '-1/10', '1/10']))


def parent(sources, state, Z=1):
    """The complete homogeneous quarter-parent; no terms omitted."""
    d, n, m, C, A, D, alpha = sources
    u, v, W, X, Y = state
    Du, Dv = 1+u*u, 1+v*v
    return (n*Du*Dv*X*X
            +m*(v*W-X+u*v*Y-u*Z)**2
            +C*Dv*(u*W-Y)**2
            +A*Du*Dv*Y*Y
            +(D+d/2)*(u*W+u*v*X-Y-v*Z)**2
            +(alpha+d/2)*Du*(v*Y-Z)**2
            +(d/2)*Dv*(X+u*Z)**2
            -(d/2)*Du*Dv*Z*Z)


def det(a):
    a = [list(row) for row in a]
    out = F(1)
    for i in range(len(a)):
        k = next((k for k in range(i, len(a)) if a[k][i]), None)
        if k is None:
            return F(0)
        if k != i:
            a[i], a[k] = a[k], a[i]
            out = -out
        p = a[i][i]
        out *= p
        for k in range(i+1, len(a)):
            ratio = a[k][i]/p
            for j in range(i+1, len(a)):
                a[k][j] -= ratio*a[i][j]
    return out


def ldl(a):
    n = len(a)
    L = [[F(i == j) for j in range(n)] for i in range(n)]
    D = []
    for j in range(n):
        dj = a[j][j] - sum((L[j][k]**2*D[k] for k in range(j)), F(0))
        assert dj > 0
        D.append(dj)
        for i in range(j+1, n):
            L[i][j] = (a[i][j]-sum((L[i][k]*L[j][k]*D[k] for k in range(j)), F(0)))/dj
    return L, D


def lower_inverse(L):
    n = len(L)
    K = [[F(i == j) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            K[i][j] = -sum((L[i][k]*K[k][j] for k in range(j, i)), F(0))
    return K


def embed(p, n, start):
    return Poly(n, {(0,)*start+e+(0,)*(n-start-p.n): c for e, c in p.d.items()})


def matrix_json(a):
    return [[str(v) for v in row] for row in a]


def write_json(name, data):
    path = ROOT/name
    path.write_text(json.dumps(data, sort_keys=True, separators=(',', ':'))+'\n', encoding='utf-8')
    return {'file': name, 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def generate():
    variables = [Poly.var(16, i) for i in range(16)]
    sources = [F(1)] + [CSTAR[i]+variables[i] for i in range(6)]
    state = [THETASTAR[i]+variables[6+i]+variables[11+i] for i in range(5)]
    Pdelta = parent(sources, state)
    contact = [Poly(11) for _ in range(6)]
    Hdelta = [[Poly(16) for _ in range(5)] for _ in range(5)]
    for e, a in sorted(Pdelta.d.items()):
        d = e[11:]
        degree = sum(d)
        if degree == 0:
            contact[0] += Poly(11, {e[:11]: a})
        elif degree == 1:
            i = d.index(1)
            contact[i+1] += Poly(11, {e[:11]: a})
        else:
            square = next((i for i, k in enumerate(d) if k >= 2), None)
            r = list(e)
            if square is not None:
                r[11+square] -= 2
                Hdelta[square][square] += Poly(16, {tuple(r): a})
            else:
                i, j = [i for i, k in enumerate(d) if k][:2]
                r[11+i] -= 1
                r[11+j] -= 1
                term = Poly(16, {tuple(r): a/2})
                Hdelta[i][j] += term
                Hdelta[j][i] += term
    assert all(not p.constant() for p in contact)
    jac = []
    for p in contact:
        row = []
        for i in range(6):
            e = [0]*11
            e[i] = 1
            row.append(p.d.get(tuple(e), F(0)))
        jac.append(row)
    minor = det(jac)
    assert minor == F(114048, 78125), minor
    H0 = [[Hdelta[i][j].constant() for j in range(5)] for i in range(5)]
    L, pivots = ldl(H0)
    assert pivots == list(map(F, ['10201/3600', '20002/91809', '27775/10001', '891/202', '1067/162']))
    K = lower_inverse(L)
    for i in range(5):
        for j in range(5):
            assert H0[i][j] == sum((L[i][k]*pivots[k]*L[j][k] for k in range(5)), F(0))
    eps = min(pivots)/20
    assert eps == F(10001, 918090)
    delta_sub = variables[:11]+[variables[11+i]-variables[6+i] for i in range(5)]
    H = [[None]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i, 5):
            H[i][j] = H[j][i] = Hdelta[i][j].substitute(delta_sub)
    Hhat = [[sum((K[i][a]*K[j][b]*H[a][b]
                   for a in range(5) for b in range(5) if K[i][a] and K[j][b]), Poly(16))
             for j in range(5)] for i in range(5)]
    linear = [tuple(L[j][i] for j in range(5)) for i in range(5)]
    frames = [(f'f{i}', linear[i], Hhat[i][i]-8*eps) for i in range(5)]
    for i in range(5):
        for j in range(i+1, 5):
            frames.append((f'f{i}+f{j}', tuple(a+b for a, b in zip(linear[i], linear[j])), eps+Hhat[i][j]/2))
            frames.append((f'f{i}-f{j}', tuple(a-b for a, b in zip(linear[i], linear[j])), eps-Hhat[i][j]/2))
    assert len(frames) == 25
    Bmax = max(g.l1_nonconstant() for _, _, g in frames)
    rho = min(F(1, 2), eps/(8*(1+Bmax)))
    expected_rho = F(510101005, 350806159660744)
    print('computed radius', rho, 'matches prior specification', rho == expected_rho, flush=True)
    assert all(rho < abs(c)/2 for c in CSTAR+THETASTAR)
    lower_bound = 5*eps/16
    frame_records = []
    blocks = []
    products = []
    product_monomials = set()
    verification = {'contact_equations_vanish_at_anchor': True,
                    'contact_minor': str(minor), 'LDL_positive': True,
                    'radius_matches_prior_specification': rho == expected_rho}
    exact_remainder = Poly(16)
    separated_remainder = Poly(16)
    intervals = []
    for frame_id, (name, coeffs, gamma) in enumerate(frames):
        gamma0 = gamma.constant()
        G, T = Poly(11), Poly(5)
        mixed = {}
        for e, c in gamma.d.items():
            a, b = e[:11], e[11:]
            if any(a) and any(b):
                mixed[b] = mixed.get(b, Poly(11)) + Poly(11, {a: c})
            elif any(a):
                G += Poly(11, {a: c})
            elif any(b):
                T += Poly(5, {b: c})
        sigma = {b: sum((abs(c) for c in a.d.values()), F(0)) for b, a in mixed.items()}
        loss = sum((sigma[b]*rho**(sum(b)+1) for b in mixed), F(0))
        p = Poly.const(11, gamma0/2-loss/2)+G
        q = Poly.const(5, gamma0/2-loss/2)+T
        mixed_records = []
        reconstructed = Poly(16)
        weights = []
        for b, a in sorted(mixed.items()):
            shift_s, shift_t = sigma[b]*rho, rho**sum(b)
            p -= a*shift_t
            q -= Poly(5, {b: shift_s})
            U = a+shift_s
            V = Poly(5, {b: F(1)})+shift_t
            mixed_records.append({'source_polynomial': U.records(),
                                  'source_constant_free_polynomial': a.records(),
                                  'source_coefficient_l1': str(sigma[b]),
                                  'source_shift': str(shift_s),
                                  'testing_monomial': list(b), 'testing_shift': str(shift_t)})
            reconstructed += embed(U, 16, 0)*embed(V, 16, 11)
            weights.append(('mixed', len(mixed_records)-1, V))
        reconstructed += embed(p, 16, 0)+embed(q, 16, 11)
        assert reconstructed == gamma
        for label, r in [('source', p), ('testing', q)]:
            bound = r.constant()-sum((abs(a)*rho**sum(e) for e, a in r.d.items() if any(e)), F(0))
            assert bound >= lower_bound, (frame_id, label, bound, lower_bound)
            intervals.append({'frame': frame_id, 'remainder': label, 'exact_coefficient_bound': str(bound)})
        deltaell = sum((a*(variables[11+i]-variables[6+i]) for i, a in enumerate(coeffs)), Poly(16))
        exact_remainder += gamma*deltaell**2
        separated_remainder += reconstructed*deltaell**2
        record = {'id': frame_id, 'name': name, 'linear_testing_function': list(map(str, coeffs)),
                  'gamma_polynomial': [[list(e[:11]), list(e[11:]), str(c)] for e, c in sorted(gamma.d.items())],
                  'gamma_residue': str(gamma0), 'nonconstant_coefficient_l1': str(gamma.l1_nonconstant()),
                  'positive_source_remainder': p.records(), 'positive_testing_remainder': q.records(),
                  'mixed_positive_products': mixed_records, 'PSD2_block_count': 2+len(mixed)}
        frame_records.append(record)
        t = [Poly.var(5, i) for i in range(5)]
        ell = sum((a*t[i] for i, a in enumerate(coeffs)), Poly(5))
        weights = [('source_remainder', None, Poly.const(5, 1)), ('testing_remainder', None, q)]+weights
        for role, mixed_id, V in weights:
            block_id = len(blocks)
            blocks.append({'id': block_id, 'frame': frame_id, 'role': role,
                           'mixed_index': mixed_id, 'testing_weight': V.records(),
                           'feature_indices': [2*block_id, 2*block_id+1]})
            for (i, j), product in zip([(0, 0), (0, 1), (1, 1)], [V, V*ell, V*ell**2]):
                for exponent, coefficient in sorted(product.d.items()):
                    products.append((block_id, i, j, exponent, coefficient))
                    product_monomials.add(exponent)
    Ptest = Pdelta.substitute(delta_sub)
    contact0 = embed(contact[0], 16, 0)
    contact1 = sum((embed(contact[i+1], 16, 0)*(variables[11+i]-variables[6+i]) for i in range(5)), Poly(16))
    assert exact_remainder == Ptest-contact0-contact1
    assert exact_remainder == separated_remainder
    monomials = sorted(product_monomials)
    monomial_index = {e: i for i, e in enumerate(monomials)}
    product_records = [[b, i, j, monomial_index[e], str(a)] for b, i, j, e, a in products]
    t = [Poly.var(5, i) for i in range(5)]
    test_state = [THETASTAR[i]+t[i] for i in range(5)]
    targets = []
    targets_full = []
    for source_id, source_name in enumerate(SOURCE_ORDER):
        unit = [F(int(source_id == j)) for j in range(7)]
        coefficient = parent(unit, test_state)
        targets_full.append({'source': source_name, 'polynomial': coefficient.records()})
        for e, c in sorted(coefficient.d.items()):
            assert e in monomial_index, (source_name, e)
            targets.append([monomial_index[e], source_id, str(c)])
    targets.sort(key=lambda r: (r[0], r[1]))
    verification.update({'frame_identity_exact': True, 'separation_identity_exact': True,
                         'all_50_positive_remainders_certified': True,
                         'all_seven_source_targets_in_product_space': True,
                         'Taylor_remainder_nonzero_monomials': len(exact_remainder.d),
                         'positive_remainder_interval_bounds': intervals})
    scope = ('Exact finite PSD(2) certificate data on the declared regular-C2 normalized '
             'contact graph and independent testing box only. No global finite cover, '
             'no global quantum-set pencil, and no global E1-E3 completion are asserted.')
    output = {'schema': 'q222-regular-c2-local-1', 'scope': scope,
              'source_order': SOURCE_ORDER, 'source_anchor': ['1']+list(map(str, CSTAR)),
              'normalizations': {'d': '1', 'Z': '1'},
              'contact_order': CONTACT_ORDER, 'contact_anchor': list(map(str, THETASTAR)),
              'source_center_variables': ['s_'+v for v in SOURCE_ORDER[1:]]+['h_'+v for v in CONTACT_ORDER],
              'testing_center_variables': ['t_'+v for v in CONTACT_ORDER],
              'contact_equations': [{'name': name, 'polynomial': p.records()} for name, p in zip(['P']+['P_'+v for v in CONTACT_ORDER], contact)],
              'contact_jacobian_source_columns': SOURCE_ORDER[1:],
              'contact_jacobian_at_anchor': matrix_json(jac), 'contact_minor_at_anchor': str(minor),
              'Taylor_allocation': 'First squared difference; otherwise first two occupied differences, half in each symmetric entry.',
              'H_delta_upper': [{'i': i, 'j': j, 'polynomial': Hdelta[i][j].records()} for i in range(5) for j in range(i, 5)],
              'H_centered_upper': [{'i': i, 'j': j, 'polynomial': H[i][j].records()} for i in range(5) for j in range(i, 5)],
              'H0': matrix_json(H0), 'LDL_L': matrix_json(L), 'LDL_D': list(map(str, pivots)),
              'LDL_L_inverse': matrix_json(K), 'epsilon': str(eps),
              'maximum_nonconstant_coefficient_l1': str(Bmax), 'radius': str(rho),
              'uniform_positive_remainder_lower_bound': str(lower_bound),
              'frames': frame_records, 'blocks': blocks,
              'fixed_PSD2_block_count': len(blocks), 'fixed_feature_count_with_repetitions': 2*len(blocks),
              'testing_product_monomials': [list(e) for e in monomials],
              'testing_product_ambient_dimension': len(monomials),
              'maximum_product_degree': max(map(sum, monomials))}
    product_output = {'schema': 'q222-regular-c2-products-1', 'scope': scope,
                      'source_order': SOURCE_ORDER, 'monomials': [list(e) for e in monomials],
                      'block_count': len(blocks), 'gram_scalar_variable_count': 3*len(blocks),
                      'coefficient_equation_count': len(monomials),
                      'symmetric_upper_product_records': product_records,
                      'off_diagonal_pairing_multiplier': 2,
                      'seven_source_target_records': targets,
                      'seven_source_target_polynomials': targets_full}
    feature_records = []
    source_gram_records = []
    for block in blocks:
        frame = frame_records[block['frame']]
        ell = Poly(5, {tuple(int(i == j) for j in range(5)): F(c)
                       for i, c in enumerate(frame['linear_testing_function']) if F(c)})
        ell_source = embed(ell, 11, 6)
        if block['role'] == 'source_remainder':
            Urecords = frame['positive_source_remainder']
        elif block['role'] == 'testing_remainder':
            Urecords = Poly.const(11, 1).records()
        else:
            Urecords = frame['mixed_positive_products'][block['mixed_index']]['source_polynomial']
        source_gram_records.append({'block': block['id'], 'nonnegative_source_weight': Urecords,
                                    'rank_one_vector': [(-ell_source).records(), Poly.const(11, 1).records()],
                                    'matrix_convention': 'U*v*v^T'})
        for component, poly in enumerate([Poly.const(5, 1), ell]):
            feature_records.append({'feature': 2*block['id']+component, 'block': block['id'],
                                    'component': component,
                                    'nonnegative_testing_radicand': block['testing_weight'],
                                    'polynomial_multiplier': poly.records(),
                                    'evaluation_convention': 'sqrt(radicand)*multiplier'})
    feature_output = {'schema': 'q222-regular-c2-features-1', 'scope': scope,
                      'source_center_variables': output['source_center_variables'],
                      'testing_center_variables': output['testing_center_variables'],
                      'features_with_repetitions': feature_records,
                      'source_gram_factorizations': source_gram_records,
                      'feature_count': len(feature_records), 'PSD2_block_count': len(blocks)}
    results = [write_json('regular_c2_frames.json', output), write_json('regular_c2_product_map.json', product_output),
               write_json('regular_c2_feature_list.json', feature_output)]
    counts = {'frames': len(frames), 'mixed_products': sum(len(f['mixed_positive_products']) for f in frame_records),
              'PSD2_blocks': len(blocks), 'listed_features': 2*len(blocks), 'ambient_monomials': len(monomials),
              'product_records': len(product_records), 'source_target_records': len(targets),
              'maximum_product_degree': max(map(sum, monomials))}
    receipt = {'status': 'PASS', 'kind': 'Fresh deterministic builder verification; not independent review',
               'scope': scope, 'counts': counts, 'verification': verification, 'generated_files': results}
    write_json('regular_c2_build_verification.json', receipt)
    print(json.dumps({'status': 'PASS', 'counts': counts, 'files': results}, indent=2), flush=True)


if __name__ == '__main__':
    generate()
