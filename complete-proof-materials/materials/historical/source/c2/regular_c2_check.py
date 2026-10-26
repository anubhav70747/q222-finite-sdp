#!/usr/bin/env python3
"""Independent exact dictionary-arithmetic audit of regular-C2 local data.

This file imports neither the builder nor any nonstandard package. It reconstructs
the complete parent, all six contact equations, every coefficient identity, every
feature, all product records, and all seven target columns. Python 3.10+.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import comb
import hashlib
import json

HERE = Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit('Run without -O: exact assertions are required.')
checks = {}


def require(condition, group, message):
    if not condition:
        raise AssertionError(group+': '+message)
    checks[group] = checks.get(group, 0)+1


def constant(n, c=1):
    c = Q(c)
    return {(0,)*n: c} if c else {}


def variable(n, i):
    e = [0]*n
    e[i] = 1
    return {tuple(e): Q(1)}


def plus(*polys):
    out = {}
    for p in polys:
        for e, c in p.items():
            out[e] = out.get(e, Q(0))+c
            if not out[e]:
                del out[e]
    return out


def scale(p, c):
    return {e: a*c for e, a in p.items() if a*c}


def times(*polys):
    if not polys:
        raise ValueError('times requires a polynomial')
    out = polys[0]
    for p in polys[1:]:
        d = {}
        for e, c in out.items():
            for f, a in p.items():
                g = tuple(x+y for x, y in zip(e, f))
                d[g] = d.get(g, Q(0))+c*a
        out = {e: c for e, c in d.items() if c}
    return out


def square(p):
    return times(p, p)


def lift(p, n, start=0):
    if not p:
        return {}
    width = len(next(iter(p)))
    return {(0,)*start+e+(0,)*(n-start-width): c for e, c in p.items()}


def parent(columns, z, n):
    d, nn, m, C, A, D, alpha = columns
    u, v, W, X, Y = z
    one = constant(n)
    Du, Dv = plus(one, square(u)), plus(one, square(v))
    return plus(
        times(nn, Du, Dv, square(X)),
        times(m, square(plus(times(v, W), scale(X, -1), times(u, v, Y), scale(u, -1)))),
        times(C, Dv, square(plus(times(u, W), scale(Y, -1)))),
        times(A, Du, Dv, square(Y)),
        times(plus(D, scale(d, Q(1, 2))), square(plus(times(u, W), times(u, v, X), scale(Y, -1), scale(v, -1)))),
        times(plus(alpha, scale(d, Q(1, 2))), Du, square(plus(times(v, Y), scale(one, -1)))),
        times(scale(d, Q(1, 2)), Dv, square(plus(X, u))),
        scale(times(d, Du, Dv), Q(-1, 2)))


def derivative(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = c*e[i]
    return out


def contact(p):
    out = {}
    for e, c in p.items():
        f = e[:6]+tuple(e[6+i]+e[11+i] for i in range(5))
        out[f] = out.get(f, Q(0))+c
    return {e: c for e, c in out.items() if c}


def replace_testing(p, contact_sign):
    """Replace the final variables by t_i + contact_sign*h_i."""
    out = {}
    for e, c in p.items():
        partial = {e[:11]+(0,)*5: c}
        for i, degree in enumerate(e[11:]):
            if not degree:
                continue
            d = {}
            for f, a in partial.items():
                for k in range(degree+1):
                    g = list(f)
                    g[11+i] += k
                    g[6+i] += degree-k
                    key = tuple(g)
                    d[key] = d.get(key, Q(0))+a*comb(degree, k)*contact_sign**(degree-k)
            partial = d
        out = plus(out, partial)
    return out


def decode(records, n):
    out = {}
    for exponent, c in records:
        e = tuple(exponent)
        require(len(e) == n and all(isinstance(x, int) and x >= 0 for x in e), 'encoding', 'bad exponent')
        a = Q(c)
        require(a != 0 and str(a) == c and e not in out, 'encoding', 'noncanonical/duplicate coefficient')
        out[e] = a
    return out


def determinant(a):
    a = [list(row) for row in a]
    out = Q(1)
    for i in range(len(a)):
        pivot = next((k for k in range(i, len(a)) if a[k][i]), None)
        if pivot is None:
            return Q(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            out = -out
        v = a[i][i]
        out *= v
        for j in range(i+1, len(a)):
            ratio = a[j][i]/v
            a[j] = [x-ratio*y for x, y in zip(a[j], a[i])]
    return out


def inverse(a):
    n = len(a)
    w = [list(a[i])+[Q(i == j) for j in range(n)] for i in range(n)]
    for i in range(n):
        k = next(k for k in range(i, n) if w[k][i])
        w[i], w[k] = w[k], w[i]
        pivot = w[i][i]
        w[i] = [x/pivot for x in w[i]]
        for j in range(n):
            if j != i:
                factor = w[j][i]
                w[j] = [x-factor*y for x, y in zip(w[j], w[i])]
    return [row[n:] for row in w]


def norm(p, n):
    return sum((abs(c) for e, c in p.items() if e != (0,)*n), Q(0))


def interval_lower(p, n, radius):
    return p.get((0,)*n, Q(0))-sum((abs(c)*radius**sum(e) for e, c in p.items() if any(e)), Q(0))


def audit():
    paths = [HERE/'regular_c2_frames.json', HERE/'regular_c2_product_map.json', HERE/'regular_c2_feature_list.json']
    frame_data, maps, feature_data = [json.loads(p.read_text()) for p in paths]
    c0 = tuple(map(Q, ['7/8', '9/8', '1/18', '11/8', '9/8', '1/18']))
    theta0 = tuple(map(Q, ['1', '1', '1', '-1/10', '1/10']))
    require(frame_data['source_order'] == ['d', 'n', 'm', 'C', 'A', 'D', 'alpha'], 'normalization', 'source order')
    require(frame_data['normalizations'] == {'d': '1', 'Z': '1'}, 'normalization', 'normalization')
    require(list(map(Q, frame_data['source_anchor'])) == [Q(1)]+list(c0), 'normalization', 'source anchor')
    require(tuple(map(Q, frame_data['contact_anchor'])) == theta0, 'normalization', 'contact anchor')
    v = [variable(16, i) for i in range(16)]
    columns = [constant(16)]+[plus(constant(16, c), v[i]) for i, c in enumerate(c0)]
    z = [plus(constant(16, c), v[11+i]) for i, c in enumerate(theta0)]
    full = parent(columns, z, 16)
    contact_expected = [contact(full)]+[contact(derivative(full, 11+i)) for i in range(5)]
    contact_actual = [decode(item['polynomial'], 11) for item in frame_data['contact_equations']]
    require(len(contact_actual) == 6, 'contact', 'equation count')
    for i, (actual, expected) in enumerate(zip(contact_actual, contact_expected)):
        require(actual == expected, 'contact', 'full contact equation '+str(i))
        require(not actual.get((0,)*11, 0), 'contact', 'anchor contact equation '+str(i))
    J = []
    for p in contact_expected:
        row = []
        for i in range(6):
            e = tuple(int(i == j) for j in range(11))
            row.append(p.get(e, Q(0)))
        J.append(row)
    require(J == [[Q(c) for c in row] for row in frame_data['contact_jacobian_at_anchor']], 'contact', 'source Jacobian')
    require(determinant(J) == Q(114048, 78125) == Q(frame_data['contact_minor_at_anchor']), 'contact', 'nonzero minor')
    Pdelta = replace_testing(full, 1)
    Hdelta_expected = [[{} for _ in range(5)] for _ in range(5)]
    for e, c in Pdelta.items():
        if sum(e[11:]) < 2:
            continue
        i = next((i for i in range(5) if e[11+i] >= 2), None)
        r = list(e)
        if i is not None:
            r[11+i] -= 2
            Hdelta_expected[i][i] = plus(Hdelta_expected[i][i], {tuple(r): c})
        else:
            i, j = [i for i in range(5) if e[11+i]][:2]
            r[11+i] -= 1
            r[11+j] -= 1
            term = {tuple(r): c/2}
            Hdelta_expected[i][j] = plus(Hdelta_expected[i][j], term)
            Hdelta_expected[j][i] = plus(Hdelta_expected[j][i], term)
    Hdelta = [[{} for _ in range(5)] for _ in range(5)]
    H = [[{} for _ in range(5)] for _ in range(5)]
    for key, matrix, expected_matrix in [('H_delta_upper', Hdelta, Hdelta_expected), ('H_centered_upper', H, None)]:
        require(len(frame_data[key]) == 15, 'Taylor', 'upper triangle size')
        pairs = set()
        for entry in frame_data[key]:
            i, j = entry['i'], entry['j']
            require(0 <= i <= j < 5 and (i, j) not in pairs, 'Taylor', 'upper triangle indexing')
            pairs.add((i, j))
            p = decode(entry['polynomial'], 16)
            expected = expected_matrix[i][j] if expected_matrix is not None else replace_testing(Hdelta[i][j], -1)
            require(p == expected, 'Taylor', key+' coefficient identity')
            matrix[i][j] = matrix[j][i] = p
    H0_expected = []
    for i in range(5):
        row = []
        for j in range(5):
            entry = derivative(derivative(full, 11+i), 11+j).get((0,)*16, Q(0))/2
            row.append(entry)
        H0_expected.append(row)
    H0 = [[Q(c) for c in row] for row in frame_data['H0']]
    L = [[Q(c) for c in row] for row in frame_data['LDL_L']]
    D = list(map(Q, frame_data['LDL_D']))
    require(H0 == H0_expected, 'LDL', 'half Hessian')
    require(all(d > 0 for d in D), 'LDL', 'positive pivots')
    require(all(L[i][i] == 1 and all(L[i][j] == 0 for j in range(i+1, 5)) for i in range(5)), 'LDL', 'unit lower triangular')
    require([[sum(L[i][k]*D[k]*L[j][k] for k in range(5)) for j in range(5)] for i in range(5)] == H0, 'LDL', 'LDL product')
    Linv = inverse(L)
    require(Linv == [[Q(c) for c in row] for row in frame_data['LDL_L_inverse']], 'LDL', 'inverse')
    eps = min(D)/20
    require(eps == Q(10001, 918090) == Q(frame_data['epsilon']), 'bounds', 'epsilon')
    Hhat = [[plus(*(scale(H[a][b], Linv[i][a]*Linv[j][b]) for a in range(5) for b in range(5))) for j in range(5)] for i in range(5)]
    coeff = [tuple(L[j][i] for j in range(5)) for i in range(5)]
    expected_frames = [(coeff[i], plus(Hhat[i][i], constant(16, -8*eps))) for i in range(5)]
    for i in range(5):
        for j in range(i+1, 5):
            for sign in [1, -1]:
                expected_frames.append((tuple(a+sign*b for a, b in zip(coeff[i], coeff[j])), plus(constant(16, eps), scale(Hhat[i][j], Q(sign, 2)))))
    Bmax = max(norm(gamma, 16) for ell, gamma in expected_frames)
    radius = min(Q(1, 2), eps/(8*(1+Bmax)))
    lower = 5*eps/16
    require(Bmax == Q(frame_data['maximum_nonconstant_coefficient_l1']), 'bounds', 'maximum coefficient norm')
    require(radius == Q(frame_data['radius']) == Q(510101005, 350806159660744), 'bounds', 'radius')
    require(lower == Q(frame_data['uniform_positive_remainder_lower_bound']), 'bounds', 'uniform positive lower bound')
    require(all(radius < abs(c)/2 for c in c0+theta0), 'bounds', 'source/contact anchor signs preserved')
    Jinv = inverse(J)
    Jpolys = [[derivative(p, j) for j in range(6)] for p in contact_expected]
    Jrelative = [[plus(*(scale(Jpolys[k][j], Jinv[i][k]) for k in range(6)),
                      constant(11, -int(i == j))) for j in range(6)] for i in range(6)]
    graph_norm = max(sum(norm(p, 11) for p in row) for row in Jrelative)
    graph_bound = graph_norm*radius
    require(all(not p.get((0,)*11, 0) for row in Jrelative for p in row), 'contact', 'relative Jacobian vanishes at anchor')
    require(graph_norm == Q(114199, 160) and graph_bound < 1, 'contact', 'contact Jacobian invertible on complete contact box')
    require(len(frame_data['frames']) == 25, 'frame', 'frame count')
    reconstructed = {}
    expected_blocks = []
    all_products = []
    all_monomials = set()
    expected_feature_data = []
    expected_grams = []
    mixed_count = 0
    for frame_id, (frame, (ell_coeff, gamma_expected)) in enumerate(zip(frame_data['frames'], expected_frames)):
        require(frame['id'] == frame_id, 'frame', 'frame numbering')
        require(tuple(map(Q, frame['linear_testing_function'])) == ell_coeff, 'frame', 'linear frame coefficients')
        gamma = decode([[a+b, c] for a, b, c in frame['gamma_polynomial']], 16)
        require(gamma == gamma_expected, 'frame', 'frame polynomial')
        gamma0 = gamma.get((0,)*16, Q(0))
        require(gamma0 == Q(frame['gamma_residue']) and gamma0 >= eps, 'frame', 'positive residue')
        require(norm(gamma, 16) == Q(frame['nonconstant_coefficient_l1']), 'frame', 'coefficient norm')
        groups = {}
        pure_s, pure_t = {}, {}
        for e, c in gamma.items():
            a, b = e[:11], e[11:]
            if any(a) and any(b):
                groups[b] = plus(groups.get(b, {}), {a: c})
            elif any(a):
                pure_s[a] = c
            elif any(b):
                pure_t[b] = c
        p = decode(frame['positive_source_remainder'], 11)
        q = decode(frame['positive_testing_remainder'], 5)
        require(interval_lower(p, 11, radius) >= lower, 'positivity', 'source remainder bound')
        require(interval_lower(q, 5, radius) >= lower, 'positivity', 'testing remainder bound')
        require(len(frame['mixed_positive_products']) == len(groups), 'separation', 'mixed group count')
        mixed_count += len(groups)
        separated = plus(lift(p, 16), lift(q, 16, 11))
        weights = [('source_remainder', None, p, constant(5)), ('testing_remainder', None, constant(11), q)]
        B = Q(0)
        ps = pure_s
        qt = pure_t
        for mixed_id, (entry, (beta, A)) in enumerate(zip(frame['mixed_positive_products'], sorted(groups.items()))):
            sigma = sum(abs(c) for c in A.values())
            source_shift, test_shift = sigma*radius, radius**sum(beta)
            U = plus(A, constant(11, source_shift))
            V = plus({beta: Q(1)}, constant(5, test_shift))
            require(decode(entry['source_constant_free_polynomial'], 11) == A, 'separation', 'mixed source coefficient')
            require(decode(entry['source_polynomial'], 11) == U, 'separation', 'shifted source factor')
            require(Q(entry['source_coefficient_l1']) == sigma and Q(entry['source_shift']) == source_shift, 'separation', 'source shift')
            require(tuple(entry['testing_monomial']) == beta and Q(entry['testing_shift']) == test_shift, 'separation', 'testing shift')
            require(interval_lower(U, 11, radius) >= 0 and interval_lower(V, 5, radius) == 0, 'positivity', 'mixed factors nonnegative')
            B += sigma*radius**(sum(beta)+1)
            ps = plus(ps, scale(A, -test_shift))
            qt = plus(qt, {beta: -source_shift})
            separated = plus(separated, times(lift(U, 16), lift(V, 16, 11)))
            weights.append(('mixed', mixed_id, U, V))
        ps = plus(ps, constant(11, gamma0/2-B/2))
        qt = plus(qt, constant(5, gamma0/2-B/2))
        require(p == ps and q == qt, 'separation', 'exact specified remainders')
        require(gamma == separated, 'separation', 'complete separated identity')
        ell = plus(*(scale(variable(5, i), c) for i, c in enumerate(ell_coeff)))
        ell_s = lift(ell, 11, 6)
        deltaell = plus(lift(ell, 16, 11), scale(lift(ell, 16, 6), -1))
        reconstructed = plus(reconstructed, times(gamma, square(deltaell)))
        require(frame['PSD2_block_count'] == len(weights), 'features', 'frame block count')
        for role, mixed_id, U, V in weights:
            block_id = len(expected_blocks)
            expected_blocks.append((frame_id, role, mixed_id, V))
            expected_grams.append((U, scale(ell_s, -1), constant(11)))
            expected_feature_data += [(block_id, 0, V, constant(5)), (block_id, 1, V, ell)]
            for (i, j), product in zip([(0, 0), (0, 1), (1, 1)], [V, times(V, ell), times(V, square(ell))]):
                for e, c in sorted(product.items()):
                    all_products.append((block_id, i, j, e, c))
                    all_monomials.add(e)
    remainder = plus(full, scale(lift(contact_expected[0], 16), -1),
                     *(scale(times(lift(contact_expected[i+1], 16), plus(v[11+i], scale(v[6+i], -1))), -1) for i in range(5)))
    require(reconstructed == remainder, 'assembly', 'full parent minus exact contact constant/gradient')
    require(len(remainder) == 1193, 'assembly', 'Taylor remainder support count')
    require(len(expected_blocks) == frame_data['fixed_PSD2_block_count'] == 785, 'features', 'block total')
    require(2*len(expected_blocks) == frame_data['fixed_feature_count_with_repetitions'] == 1570, 'features', 'feature total')
    require(len(frame_data['blocks']) == len(expected_blocks), 'features', 'block record count')
    for i, (record, expected) in enumerate(zip(frame_data['blocks'], expected_blocks)):
        frame_id, role, mixed_id, V = expected
        require((record['id'], record['frame'], record['role'], record['mixed_index'], record['feature_indices']) == (i, frame_id, role, mixed_id, [2*i, 2*i+1]), 'features', 'block metadata')
        require(decode(record['testing_weight'], 5) == V, 'features', 'block radicand')
    monomials = sorted(all_monomials)
    index = {e: i for i, e in enumerate(monomials)}
    require(len(monomials) == 237 and max(map(sum, monomials)) == 6, 'products', 'ambient product space')
    require(monomials == list(map(tuple, frame_data['testing_product_monomials'])) == list(map(tuple, maps['monomials'])), 'products', 'monomial indexing')
    expected_records = [[b, i, j, index[e], str(c)] for b, i, j, e, c in all_products]
    require(expected_records == maps['symmetric_upper_product_records'], 'products', 'all product coefficients')
    require(len(expected_records) == 27747, 'products', 'product record count')
    require(maps['off_diagonal_pairing_multiplier'] == 2, 'products', 'symmetric trace pairing')
    require(maps['gram_scalar_variable_count'] == 2355 and maps['coefficient_equation_count'] == 237, 'products', 'linear system dimensions')
    test = [plus(constant(5, c), variable(5, i)) for i, c in enumerate(theta0)]
    target_records = []
    for source_id in range(7):
        cols = [constant(5, int(j == source_id)) for j in range(7)]
        coeffs = parent(cols, test, 5)
        require(decode(maps['seven_source_target_polynomials'][source_id]['polynomial'], 5) == coeffs, 'targets', 'complete parent source column')
        require(all(e in index for e in coeffs), 'targets', 'source column supported')
        target_records += [[index[e], source_id, str(c)] for e, c in coeffs.items()]
    target_records.sort(key=lambda r: (r[0], r[1]))
    require(target_records == maps['seven_source_target_records'] and len(target_records) == 295, 'targets', 'all seven source records')
    require(len(feature_data['features_with_repetitions']) == len(expected_feature_data), 'explicit_features', 'feature record count')
    require(len(feature_data['source_gram_factorizations']) == len(expected_grams), 'explicit_features', 'Gram record count')
    for i, (record, expected) in enumerate(zip(feature_data['features_with_repetitions'], expected_feature_data)):
        b, component, V, multiplier = expected
        require((record['feature'], record['block'], record['component']) == (i, b, component), 'explicit_features', 'feature index')
        require(decode(record['nonnegative_testing_radicand'], 5) == V, 'explicit_features', 'feature radicand')
        require(decode(record['polynomial_multiplier'], 5) == multiplier, 'explicit_features', 'feature multiplier')
    for i, (record, (U, x, y)) in enumerate(zip(feature_data['source_gram_factorizations'], expected_grams)):
        require(record['block'] == i, 'explicit_features', 'Gram index')
        require(decode(record['nonnegative_source_weight'], 11) == U, 'explicit_features', 'Gram weight')
        require([decode(p, 11) for p in record['rank_one_vector']] == [x, y], 'explicit_features', 'Gram vector')
    receipt = {'status': 'PASS',
               'scope': 'Regular-C2 local normalized contact graph/testing box only; global E1-E3 remain open.',
               'method': 'Independent standard-library Fraction dictionary implementation; no builder import.',
               'contact_graph_box_bound': {
                   'relative_Jacobian_row_coefficient_norm': str(graph_norm),
                   'relative_Jacobian_infinity_norm_upper_bound': str(graph_bound),
                   'interpretation': 'Less than one implies source Jacobian invertible on the entire declared contact box.'},
               'counts': {'frames': 25, 'mixed_products': mixed_count, 'PSD2_blocks': len(expected_blocks),
                          'listed_features': len(expected_feature_data), 'ambient_monomials': len(monomials),
                          'product_records': len(expected_records), 'source_target_records': len(target_records),
                          'Taylor_remainder_monomials': len(remainder)},
               'assertions_by_group': checks, 'total_assertions': sum(checks.values()),
               'files_checked': [{'file': p.name, 'bytes': p.stat().st_size,
                                  'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]}
    (HERE/'regular_c2_independent_verification.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    audit()
