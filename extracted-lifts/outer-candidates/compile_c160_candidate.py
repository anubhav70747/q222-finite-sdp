"""Exact compilation of the recovered C160 OUTER CANDIDATE.

This script proves the coefficient identities and evaluates actual rational
pencils. It does not prove equality with the quantum set. Python 3 + SymPy.
"""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
from fractions import Fraction as F
import sympy as s

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'C160_source_compiler_cache.json'
OUT = ROOT / 'generated_c160'
z = s.symbols('k a0 b0 a1 b1')
k, a0, b0, a1, b1 = z
EXP = list(itertools.product(range(3), repeat=5))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(name, obj):
    p = OUT / name
    p.write_text(json.dumps(obj, separators=(',', ':')) + '\n')
    return p


def rational(q):
    q = s.Rational(q)
    return int(q.p), int(q.q)


def poly(expr):
    return s.Poly(expr, *z, domain=s.QQ)


def pterms(p):
    return [[list(e), *rational(c)] for e, c in p.terms()]


def map_from_cache(d):
    refs = d['tensor_phi_basis']['class_reference_coords']
    refmap = {c: (r, s.Rational(w)) for r, (c, w) in enumerate(refs)}
    cols = {c: {r: w} for c, (r, w) in refmap.items()}
    for v in d['nullspace_basis']:
        free = v['free']
        assert sum(s.Rational(n, den) for c, n, den in v['entries'] if c == free) == 1
        col = {}
        for c, n, den in v['entries']:
            if c == free:
                continue
            assert c in refmap
            r, w = refmap[c]
            col[r] = col.get(r, 0) - s.Rational(n, den) * w
        cols[free] = {r: w for r, w in col.items() if w}
    assert set(cols) == set(range(2140))
    T = s.zeros(243, 9)
    for j, vec in enumerate(d['target_affine_section']):
        for c, n, den in vec:
            for r, w in cols[c].items():
                T[r, j] += s.Rational(n, den) * w
    return cols, T


def sparse_linear_mat_entries(d, cols):
    """Raw adjoint matrices: divide actual off-diagonal columns by two."""
    records = []
    for b in d['block_layout']:
        c = b['coord_start']
        for i in range(b['size']):
            for j in range(i, b['size']):
                for r, w in sorted(cols[c].items()):
                    records.append((b['block'], i, j, r, w / (1 if i == j else 2)))
                c += 1
        assert c == b['coord_end']
    return records


def check_and_export_features(d, cols, T):
    phi = []
    for e in EXP:
        p = s.Integer(1)
        for x, a in zip(z, e):
            p *= ((1-x)**2, x*(1-x), x*x)[a]
        phi.append(poly(p))
    P = {c: sum((w*phi[r] for r, w in col.items()), poly(0)) for c, col in cols.items()}
    gaps = [k*b0*(1-a0)-a0*(1-b0), k*a1*(1-b0)-b0*(1-a1), k*b1*(1-a1)-a1*(1-b1)]
    atoms = [poly(x) for x in z] + [poly(1-x) for x in z] + [poly(g) for g in gaps]
    names = [str(x) for x in z] + ['1-'+str(x) for x in z] + ['g1','g2','g3']
    features = []
    identities = 0
    for b in d['block_layout']:
        free = b['free']
        bits = list(itertools.product(range(2), repeat=len(free)))
        v = [poly(s.prod(z[f] if bit else 1-z[f] for f, bit in zip(free, bitvec))) for bitvec in bits]
        assert len(v) == b['size']
        w, rem = P[b['coord_start']].div(v[0]*v[0])
        assert rem.is_zero
        remaining = w
        factors = []
        for name, atom in zip(names, atoms):
            power = 0
            while True:
                quotient, rem = remaining.div(atom)
                if not rem.is_zero:
                    break
                remaining = quotient
                power += 1
            if power:
                factors.append([name, power])
        assert remaining.total_degree() == 0 and remaining.LC() > 0, (b['block'], remaining)
        c = b['coord_start']
        for i in range(b['size']):
            for j in range(i, b['size']):
                assert P[c] == (1 if i == j else 2)*w*v[i]*v[j], (b['block'], i, j)
                identities += 1
                c += 1
        features.append({'block':b['block'], 'size':b['size'], 'free_variable_indices':free,
                         'basis_bits':bits, 'weight_factors':factors,
                         'weight_constant':list(rational(remaining.LC())), 'weight_power_coefficients':pterms(w)})
    dz = lambda x: (1-x)**2+x*x
    D = (1+k*k)*s.prod(dz(x) for x in z[1:])
    ca = [(1-2*x)/dz(x) for x in (a0,a1)]
    sa = [2*x*(1-x)/dz(x) for x in (a0,a1)]
    cb = [(1-2*x)/dz(x) for x in (b0,b1)]
    sb = [2*x*(1-x)/dz(x) for x in (b0,b1)]
    q = [(1-k*k)/(1+k*k)*c for c in ca+cb]
    q += [ca[x]*cb[y]+2*k/(1+k*k)*sa[x]*sb[y] for x in range(2) for y in range(2)]
    physical = [poly(D)] + [poly(s.cancel(D*v)) for v in q]
    for j in range(9):
        target = sum((T[r,j]*phi[r] for r in range(243) if T[r,j]), poly(0))
        assert target == (1 if j == 0 else -1)*physical[j], j
    # This positive definite Gram normalizes the dual and controls all moments.
    QD = s.eye(16).row_join(s.eye(16)).col_join(s.eye(16).row_join(2*s.eye(16)))
    got = s.zeros(243,1)
    c = 0
    for i in range(32):
        for j in range(i,32):
            if QD[i,j]:
                for r,w in cols[c].items():
                    got[r] += w*QD[i,j]
            c += 1
    assert got == T[:,0]
    # Exact strict-grid bounds also verify that every literal feature weight is nonzero somewhere.
    center = [s.Rational(1,2),s.Rational(1,100),s.Rational(1,10),s.Rational(1,2),s.Rational(9,10)]
    grid = [[x+s.Rational(e,1000) for x,e in zip(center,eps)] for eps in itertools.product([-1,1],repeat=5)]
    min_gap = min(g.subs(dict(zip(z,point))) for point in grid for g in gaps)
    assert min_gap > 0
    out = {'classification':'FIXED FEATURES FOR AN OUTER CANDIDATE; GLOBAL COMPLETENESS UNPROVED',
           'variables':list(map(str,z)), 'domain':'0<=z_i<=1; g1>=0,g2>=0,g3>=0',
           'gaps':[pterms(poly(g)) for g in gaps],
           'common_denominator':pterms(poly(D)), 'denominator_bounds':['1/16','2'],
           'physical_numerators':[pterms(p) for p in physical[1:]],
           'basis_rule':'v_b(z)=sqrt(weight_b(z)/D(z)) times product of z_i or (1-z_i) over free variables, in basis_bits order',
           'blocks':features, 'total_feature_count':sum(b['size'] for b in d['block_layout']),
           'feature_product_identity_count':identities,
           'strict_grid_min_gap':list(rational(min_gap)),
           'normalization_gram':'block 0: [[I16,I16],[I16,2I16]]',
           'source_cache_sha256':sha(CACHE)}
    dump('C160_FIXED_FEATURES.json',out)
    return identities


def compile_canonical(d, cols, T):
    L = T.T.copy()
    for i in range(1,9):
        L[i,:] *= -1
    _, pivots = L.rref()
    assert len(pivots) == 9
    free = [i for i in range(243) if i not in pivots]
    C = L[:,pivots].inv()
    S = s.zeros(243,243)
    for a,p in enumerate(pivots):
        S[p,:9] = C[a,:]
    for j,r in enumerate(free):
        S[r,9+j] = 1
        vals = -C*L[:,r]
        for a,p in enumerate(pivots):
            S[p,9+j] = vals[a]
    assert L*S == s.eye(9).row_join(s.zeros(9,234))
    entries = {}
    for b,i,j,r,w in sparse_linear_mat_entries(d,cols):
        for var in range(243):
            if S[r,var]:
                key = (b,i,j,var)
                entries[key] = entries.get(key,0)+w*S[r,var]
    coo = [[*key,*rational(w)] for key,w in sorted(entries.items()) if w]
    obj = {'classification':'EXACT RATIONAL OUTER CANDIDATE FOR CANONICAL ATOMS; REVERSE INCLUSION UNPROVED',
           'coordinate_order':['a0','a1','b0','b1','E00','E01','E10','E11'],
           'matrix_order':499,'visible_count':8,'auxiliary_count':234,'block_count':160,'maximum_block_order':32,
           'blocks':[[b['block'],b['size']] for b in d['block_layout']],
           'format':'[block,row,col,variable,numerator,denominator]; upper triangle; mirror once; variable0=constant; 1..8=visible; 9..242=auxiliary',
           'entries':coo,'moment_pivot_indices':list(pivots),'moment_free_indices':free,
           'moment_reconstruction':[[r,j,*rational(S[r,j])] for r in range(243) for j in range(243) if S[r,j]],
           'source_cache_sha256':sha(CACHE)}
    dump('C160_CANONICAL_OUTER_PENCIL.coo.json',obj)
    return L, len(coo)


def compile_global(d, cols, L):
    vertices = []
    for aa0,aa1,bb0,bb1 in itertools.product([-1,1],repeat=4):
        vertices.append([aa0,aa1,bb0,bb1,aa0*bb0,aa0*bb1,aa1*bb0,aa1*bb1])
    V = s.Matrix([[1]+v for v in vertices]).T
    _, dependent = V.rref()
    assert len(dependent) == 9
    independent = [i for i in range(16) if i not in dependent]
    C = V[:,dependent].inv()
    actions = [s.Matrix(h) for h in d['bell_target_actions']]
    assert len(actions) == 128 and len({tuple(h) for h in actions}) == 128
    for h in actions:
        assert h*h.T == s.eye(8)
        assert {tuple(h.T*s.Matrix(v)) for v in vertices} == {tuple(v) for v in vertices}
    raw = sparse_linear_mat_entries(d,cols)
    entries = []
    for g in range(128):
        entries.extend([[160*g+b,i,j,9+243*g+r,*rational(w)] for b,i,j,r,w in raw])
    dep_linear = {idx:{j:C[a,j] for j in range(9) if C[a,j]} for a,idx in enumerate(dependent)}
    for g,h in enumerate(actions):
        action9 = s.diag(1,h.T)
        R = -C*action9*L
        for a,idx in enumerate(dependent):
            for r in range(243):
                if R[a,r]:
                    dep_linear[idx][9+243*g+r] = R[a,r]
    for j,idx in enumerate(independent):
        var = 9+128*243+j
        entries.append([128*160+idx,0,0,var,1,1])
        R = -C*V[:,idx]
        for a,dep in enumerate(dependent):
            if R[a]:
                dep_linear[dep][var] = R[a]
    for idx,linear in dep_linear.items():
        entries.extend([[128*160+idx,0,0,var,*rational(w)] for var,w in sorted(linear.items())])
    blocks = [[160*g+b['block'],b['size']] for g in range(128) for b in d['block_layout']]
    blocks.extend([[128*160+i,1] for i in range(16)])
    obj = {'classification':'GLOBAL QUANTUM OUTER CANDIDATE; EQUALITY WITH Q222 IS NOT PROVED',
           'coordinate_order':['a0','a1','b0','b1','E00','E01','E10','E11'],
           'matrix_order':63888,'visible_count':8,'auxiliary_count':31111,'block_count':20496,'maximum_block_order':32,
           'format':'[block,row,col,variable,numerator,denominator]; upper triangle; mirror once; variable0=constant; 1..8=visible; 9..31119=auxiliary',
           'blocks':blocks, 'entries':entries,
           'auxiliary_order':{'first_31104':'128 successive moment vectors of length243; exponent order product(range(3),repeat=5)',
                              'last_7_free_local_vertices':independent},
           'eliminated_local_vertices':list(dependent),'local_vertices':vertices,
           'elimination_inverse':[[str(x) for x in row] for row in C.tolist()],
           'bell_actions_on_primal_coordinates':[h.T.tolist() for h in actions],
           'source_cache_sha256':sha(CACHE)}
    # JSON cannot serialize SymPy Integer directly.
    obj['bell_actions_on_primal_coordinates'] = [[[int(v) for v in row] for row in h.T.tolist()] for h in actions]
    assert sum(order for _,order in blocks) == obj['matrix_order']
    assert len(blocks) == obj['block_count']
    assert 128*243+7 == obj['auxiliary_count']
    assert max(row[3] for row in entries) == 8+obj['auxiliary_count']
    dump('C160_GLOBAL_OUTER_PENCIL.coo.json',obj)
    return len(entries), list(dependent), independent


def main():
    OUT.mkdir(exist_ok=True)
    assert sha(CACHE) == 'd17f367bcee162479db4a31ac86808c04198bb15aa8b6ec2bdf354096d3ef703'
    d = json.loads(CACHE.read_text())
    cols,T = map_from_cache(d)
    feature_count = check_and_export_features(d,cols,T)
    print('PASS: all2140 literal feature product identities, nine physical targets, normalization Gram',flush=True)
    L,canonical_nnz = compile_canonical(d,cols,T)
    print('PASS: canonical moment elimination and evaluated pencil',flush=True)
    global_nnz,dep,indep = compile_global(d,cols,L)
    print('PASS: global affine elimination, Bell group, and evaluated pencil',flush=True)
    receipt={'scope':'Exact compiler verification. Quantum reverse inclusion was NOT proved.',
             'source_cache_sha256':sha(CACHE), 'status':'PASS_WITH_OPEN_REVERSE_INCLUSION',
             'groups':{'literal_feature_products':feature_count,'physical_target_polynomials':9,
                       'normalization_Gram':'PASS','canonical_moment_elimination':'PASS',
                       'Bell_actions':128,'global_affine_elimination':'PASS'},
             'canonical_entry_count':canonical_nnz,'global_entry_count':global_nnz,
             'global_counts':{'matrix_order':63888,'auxiliary_count':31111,'blocks':20496},
             'open_statement':'Every point in the global projected spectrahedron is in Q222.',
             'files':{p.name:sha(p) for p in sorted(OUT.glob('*.json')) if p.name != 'C160_EXTRACTION_VERIFICATION.json'}}
    dump('C160_EXTRACTION_VERIFICATION.json',receipt)
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
