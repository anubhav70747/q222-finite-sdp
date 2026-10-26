"""One exact genuine-source obstruction to the task3 fixed scalar completion.

Standalone SymPy checker. No solver, samples, package import or network.
Reconstructs the correct interpolation, including the Y_u exception.
"""
from pathlib import Path
import argparse
import json
import sympy as s

u, v, tau, W, X, Y, Z = s.symbols('u v tau W X Y Z', real=True)
Q = s.Rational
a, b, c, t, r = 1, 1, 1, Q(3, 4), Q(1, 4)
h = b-a*a*t
g = (1+c*c)*b-(1-r+a*a*c*c*t)
j = 1+c*c*r-t
Delta = (1+(1+c*c)*b+c*c*(1+c*c)*b*b+2*c*c*(1+c*c)*b*r
         +c**4*r*r+2*c*c*r+a*a*c*c*((1+c*c)*b+1+c*c*r*(2-r)))
Ns = (1+c*c)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*c*c)*(1-r)*t
Nn = (b+r)*(1+c*c*r)*(1+c*c*r+a*a*c*c*t)-r*Delta
gaps = [h, g, j, 1-c*c*h, b-a*a*c*c*r, Ns, Nn, Delta]
assert all(q > 0 for q in gaps) and 0 < r < 1
n = Nn/(2*r*Delta)
m = (1+a*a*c*c)*(1-r)*(1+c*c*r)*(1+c*c*r+a*a*c*c*t)/(2*g*Delta)
C = (1+c*c*r)*Ns/(2*a*a*t*Delta)
A = a*a*(1+c*c*r)*(b+r)*(1+c*c*r+a*a*c*c*t)/(2*h*Delta)
D = (1+c*c)*(b+r)*(b-a*a*c*c*r)*(1+c*c*r+a*a*c*c*t)/(2*a*a*j*Delta)
alpha = c**4*(b+r)*Ns/(2*(1-c*c*h)*Delta)
weights = [n, m, C, A, D, alpha]
assert weights == [Q(1, 8), Q(3, 8), Q(1, 6), Q(5, 8), Q(3, 8), Q(1, 6)]

rows, k2 = [], []
for f in (1, u):
    for g0 in (1, v):
        rows.append([0, f*g0, 0, 0]); k2.append(4*n)
rows.append([v, -1, u*v, -u]); k2.append(4*m)
for g0 in (1, v):
    rows.append([u*g0, 0, -g0, 0]); k2.append(4*C)
for f in (1, u):
    for g0 in (1, v):
        rows.append([0, 0, f*g0, 0]); k2.append(4*A)
rows.append([u, u*v, -1, -v]); k2.append(4*D+2)
for f in (1, u):
    rows.append([0, 0, f*v, -f]); k2.append(4*alpha+2)
for g0 in (1, v):
    rows.append([0, g0, 0, u*g0]); k2.append(Q(2))
M = s.Matrix(rows)
R = M.col_join(v*M)
kk = k2+k2
S = s.zeros(32, 29)
col = 0
for i in (1, 2, 3, 14, 15):
    S[i, col] = 1; S[0, col] = -M[i, 1]; col += 1
for i in (8, 9, 10, 12, 13):
    S[i, col] = 1; S[7, col] = -M[i, 2]; col += 1
for entries in ({6: 1, 5: -v}, {11: 1, 5: -1, 0: -u*v},
                {4: u, 5: -v, 0: u, 7: -v, 9: -u*v}):
    for i, value in entries.items():
        S[i, col] = value
    col += 1
for i in range(16):
    S[16+i, col] = 1; S[i, col] = -v; col += 1
high = s.zeros(32, 1)
for i, value in {27: v*v, 20: -u*v, 31: -u*v*v,
                 16: -u*v, 24: v, 26: u*v}.items():
    high[i] = value
assert all(s.expand(x) == 0 for x in R[:, :3].T*S)
assert all(s.expand(x) == 0 for x in R[:, :3].T*high)
contact = {u: 1, v: 1}
psi_c = s.Matrix([1, -Q(1, 4), Q(1, 4), 1])
e_norm = R.subs(contact)*psi_c
bc = s.Matrix([-kk[i]*e_norm[i] for i in range(32)])

def interpolate(bb):
    eta = [bb[i]+bb[16+i] for i in range(16)]
    indices = (1, 2, 3, 14, 15, 8, 9, 10, 12, 13, 6, 11, 4)
    # Here u_c=v_c=1. The Y_u coefficient must include eta4.
    values = [eta[i] if i != 9 else eta[9]+eta[4] for i in indices]
    return s.Matrix(values+list(bb[16:]))

assert s.simplify(S.subs(contact)*interpolate(bc)-bc) == s.zeros(32, 1)
ww = S*interpolate(bc)+tau*(high-S*interpolate(high.subs(contact)))
scalar = (-2*(1+u*u)*(1+v*v)**2
          -2*sum(R[i, 3]*ww[i] for i in range(32))
          -sum(ww[i]**2/kk[i] for i in range(32)))
line = s.factor(scalar.subs(v, 1))
expected = -(u-1)**2*(379*tau*tau-1476*tau+2217)/30
assert s.expand(line-expected) == 0
positive_bracket = 379*(tau-Q(738, 379))**2+Q(295599, 379)
assert s.expand(positive_bracket-(379*tau*tau-1476*tau+2217)) == 0
assert s.factor(line.subs({u: Q(1, 2), tau: 2})) == -Q(781, 120)
assert s.factor(s.diff(line, u, 2).subs({u: 1, tau: 2})) == -Q(781, 15)

psi = s.Matrix([W, X, Y, Z])
full = sum(kk[i]*(R[i, :]*psi)[0]**2 for i in range(32))-2*(1+u*u)*(1+v*v)**2*Z*Z
cs = {u: 1, v: 1, W: 1, X: -Q(1, 4), Y: Q(1, 4), Z: 1}
assert all(s.factor(f.subs(cs)) == 0 for f in [full]+[s.diff(full, z) for z in (W, X, Y, Z, u, v)])
fixed_state_line = s.factor(full.subs({v: 1, W: 1, X: -Q(1, 4), Y: Q(1, 4), Z: 1}))
assert fixed_state_line == Q(289, 24)*(u-1)**2
jet_norm_squared = s.factor(fixed_state_line-line.subs(tau, 2))
assert jet_norm_squared == Q(1523, 40)*(u-1)**2
point = {u: Q(1, 2), v: 1, tau: 2}
completed_sum = s.factor(sum(kk[i]*(R[i, 3]+ww[i]/kk[i])**2 for i in range(32)).subs(point))
assert completed_sum == Q(517, 40)
assert s.factor(full.subs({**point, W: 0, X: 0, Y: 0, Z: 1})) == Q(77, 12)
assert Q(77, 12)-completed_sum == -Q(781, 120)
out = {
    'status': 'EXACT OBSTRUCTION TO THE SELECTED SCALAR COMPLETION; NOT TO K364',
    'native_source_a_b_c_t_r': [str(x) for x in (a, b, c, t, r)],
    'gaps_h_g_j_E_Q_Nstar_Nn_Delta': [str(x) for x in gaps],
    'weights_n_m_C_A_D_alpha': [str(x) for x in weights],
    'contact': ['1', '1', '1', '-1/4', '1/4', '1'],
    'scalar_line_all_real_tau': str(line),
    'strictly_positive_bracket': str(positive_bracket),
    'tau2_test_u_v': ['1/2', '1'],
    'tau2_scalar_at_test': '-781/120',
    'tau2_scalar_u_second_derivative_at_contact': '-781/15',
    'tau2_highest_Bob_pivot': '6/7',
    'full_DvF_at_test_state_0001': '77/12',
    'completed_squares_at_test_state_0001': '517/40',
    'fixed_contact_state_line_DvF': str(fixed_state_line),
    'current_completed_root_jet_squared': str(jet_norm_squared),
    'genuine_source_stationarity': 'value and all six state/angular derivatives zero exactly',
    'remaining': 'Additional zero-contact syzygy directions or a different K364 Gram; H2 remains open.'
}
ap = argparse.ArgumentParser()
ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
args = ap.parse_args()
args.output.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
print(json.dumps(out, indent=2))
