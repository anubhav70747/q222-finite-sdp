#!/usr/bin/env python3
"""Cross-review corroboration of diagonal/B0 records, using exact Fractions.

The analytic review observations are stated explicitly in the receipt and are
not described as computer-proved. This script does not import the diagonal
generator. It imports only arithmetic from this folder's independent checker.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
from regular_c2_check import constant, variable, plus, times, square, scale, lift, derivative

HERE = Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit('Run without -O: exact assertions are required.')
BASE = HERE.parent
DIAG = BASE/'diagonal'


def readpoly(record):
    return {tuple(t['powers']): Q(t['coefficient']) for t in record['terms'] if Q(t['coefficient'])}


def evaluate(p, values):
    return sum((c*__import__('functools').reduce(lambda a, b: a*b,
                 (Q(v)**k for v, k in zip(values, e)), Q(1)) for e, c in p.items()), Q(0))


def run():
    b0_path = DIAG/'B0_EXACT.json'
    event_path = DIAG/'EVENT_CONTACT_EXACT.json'
    b0, event = [json.loads(p.read_text()) for p in [b0_path, event_path]]
    count = 0
    U, V = variable(6, 0), variable(6, 1)
    z = [variable(6, i+2) for i in range(4)]
    one = constant(6)
    omega = plus(times(U, z[0]), scale(z[1], -1), scale(times(U, V, z[3]), -1))
    target = [times(plus(one, square(U)), square(z[0])),
              times(plus(one, square(V)), square(omega)),
              times(plus(one, square(U)), plus(one, square(V)), square(z[1])),
              times(square(U), square(z[2]))]
    features = []
    permitted = {(0, 0), (1, 0), (0, 1), (1, 1), (1, 2)}
    for row in b0['feature_coefficients_in_state_order_z0_z1_z2_z3']:
        feature = {}
        for i, entry in enumerate(row):
            p = readpoly(entry)
            assert set(p).issubset(permitted)
            count += 1
            feature = plus(feature, times(lift(p, 6), z[i]))
        features.append(feature)
    assert len(features) == 9
    count += 1
    for i, label in enumerate(['c1', 'c2', 'c3', 'c4']):
        actual = plus(*(square(f) for f, w in zip(features, b0['weights']) if w == label))
        assert actual == target[i]
        count += 1
    expected_tests = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 4, 0], [0, 0, 0, 1]]
    for record, expected in zip(b0['necessary_tests'], expected_tests):
        values = [record['U'], record['V']]+record['state']
        assert [evaluate(p, values) for p in target] == list(map(Q, expected))
        count += 1
    T, L, A, eps0 = Q(1), Q(2, 5), Q(1, 5), Q(1)
    p, w, h = 1-A-L*T*T, 1-L*T*T, 1-L*(1+T*T)
    coefficients = [A*T*(1+T*T)*p/w, A*eps0, eps0*p, A*eps0*(1+T*T)*w*p/h]
    assert list(map(str, [p, w, h])) == b0['sample_p_w_h']
    assert list(map(str, coefficients)) == b0['sample_coefficients']
    assert L*(T*T+2)-1 == Q(1, 5)
    count += 3
    assert b0['full_B0_restoration_supplied'] is False
    assert b0['congruence_preserves_order_two_blocks_claimed'] is False
    count += 2
    b0_count = count
    u, v, a, b = [variable(4, i) for i in range(4)]
    one4 = constant(4)
    Du, Dv = plus(one4, square(u)), plus(one4, square(v))
    forms = [scale(plus(times(plus(Dv, times(square(u), square(v))), square(a)),
                        scale(times(square(u), v, a), -2), square(u), times(Dv, square(b)),
                        scale(times(u, Dv, b), 2)), 2),
             scale(times(Du, Dv, square(b)), 4),
             scale(square(plus(times(u, v, a), scale(u, -1), scale(times(Dv, b), -1))), 4),
             scale(times(Dv, square(plus(a, times(u, v, b)))), 4),
             scale(times(Du, Dv, square(a)), 4),
             scale(square(plus(a, v)), 4)]
    assert [readpoly(p) for p in event['forms']] == forms
    count += 1
    jac = [forms]+[[derivative(p, i) for p in forms] for i in range(4)]
    assert [[readpoly(p) for p in row] for row in event['jacobian_rows']] == jac
    count += 1
    x = [readpoly(p) for p in event['primitive_cofactor_components']]
    assert [max(map(sum, p)) for p in x] == [15, 13, 13, 11, 15, 17]
    assert [len(p) for p in x] == [59, 40, 43, 36, 42, 38]
    count += 2
    for row in jac:
        assert not plus(*(times(p, c) for p, c in zip(row, x)))
        count += 1
    reviewed = [b0_path, event_path, DIAG/'diagonal_proofs.tex',
                BASE/'inputs/frozen_august31/Q222_SXDEG2_CORRECTION_PATCH_AND_FREEZE_2026-08-31.md',
                BASE/'inputs/canonical_august30/Q222_SXDEG2_CONSOLIDATED_PARCHMENT_FINAL_2026-08-30.md']
    report = {
        'status': 'PASS_WITH_STATED_LOCAL_AND_CONDITIONAL_SCOPE',
        'reviewer': 'Separate C2 reconstruction agent; same model-assisted system, not external review.',
        'algebraic_corroboration': {
            'method': 'Exact Fraction polynomial dictionaries; no diagonal-generator import.',
            'B0_checks': b0_count, 'event_checks': count-b0_count,
            'B0_result': 'All nine feature functions reproduce all four complete coefficient columns; four rational necessity tests and the strict-chamber sample agree.',
            'event_result': 'All six displayed restricted forms, all thirty derivative entries, the degree/support tuples, and all five primitive-vector kernel identities agree.',
            'event_limit': 'This cross-review does not independently recompute the six determinants or certify source physicality.'},
        'analytic_review': [
            {'topic': 'Positive separated coefficient identity',
             'finding': 'Expansion and the stated strict p/q lower bounds are correct on the required product neighborhood.'},
            {'topic': 'Smooth-diagonal normal correction',
             'finding': 'A*Delta z=0 makes the normal correction exact; positive tangential H and positive normal A^T A give a positive full matrix for sufficiently large integer t by a Schur complement.'},
            {'topic': 'Constant frame change and finite frame',
             'finding': 'For H*=C^T C and f=Cz, C^-T(H+tA^T A)C^-1 has identity residue. The coefficients K_ii-2(m-1)eps and eps+/-K_ij/2 reproduce every entry.'},
            {'topic': 'Mixed denominator warning',
             'finding': 'The independent-poles argument proves infinite separation rank of1/(1-xy) on product intervals; the construction correctly excludes unsupported mixed inversions.'},
            {'topic': 'Frozen Repair5 factor typing',
             'finding': 'August30 equation12.1 and August31 Repair5 put A,T,p,w,h,eps0 in the point/coefficient factor and U,V in the support factor. The new proof retains this and labels its testing notation as a coefficient-module abstraction.'},
            {'topic': 'B0 restoration scope',
             'finding': 'The order4 diagonal orthant pencil concerns only(c1,c2,c3,c4). The20-function congruence envelope assumes a pure coefficient-factor state change. No physical eight-coordinate pencil or surviving order2 bound is inferred. The missing congruence and clearing maps are explicit.'},
            {'topic': 'Finite-cover indicator assembly',
             'finding': 'Membership-pattern cells are finite and choose a covering rectangle for every nonempty source/test cell pair. Zero extension preserves fixed feature size and yields the asserted global slack identity under the supplied finite-cover hypothesis.'},
            {'topic': 'Moment product space and coordinates',
             'finding': 'W explicitly contains1, all Phi_i, and every feature product. Full dimensionality justifies extending1,Phi_i to a basis. No unsupported product-span assumption remains.'},
            {'topic': 'Both inclusions and closure',
             'finding': 'Forward inclusion uses an actual finite convex combination of evaluations. Reverse inclusion uses separation of the assumed closed C=convPhi(T) and the exact global slack identities. No closure of a projected spectrahedron is asserted or required.'}],
        'remaining_global_inputs': [
            'An evaluated finite all-source/all-test cover with restored local slack identities.',
            'All exceptional, rank-drop, and projective boundary coverage and restoration maps.',
            'An evaluated global feature basis, coefficient tables, and affine projection in the original eight quantum coordinates.'],
        'files_reviewed': [{'file': str(p.relative_to(BASE)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in reviewed]}
    target_path = HERE/'regular_c2_diagonal_cross_review.json'
    target_path.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': report['status'], 'B0_checks': b0_count, 'event_checks': count-b0_count,
                      'receipt': target_path.name}, indent=2))


if __name__ == '__main__':
    run()
