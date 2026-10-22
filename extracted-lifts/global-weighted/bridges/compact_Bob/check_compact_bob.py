#!/usr/bin/env python3
"""Exact new identities only; does not run or regenerate the 669/377 matrices.

Requires SymPy.  Run from any directory.  The companion JSON records the
small evaluated coefficient matrices and the fixed-family prescription.
"""
from pathlib import Path
import json
import sympy as S

L, r, u, c, s = S.symbols("L r u c s")
P, X, Y, T, Z = S.symbols("P X Y T Z")
W = S.symbols("W")
C, A, a = S.symbols("C A a")
Cp, Ap, ap = S.symbols("Cprime Aprime aprime")
h = S.Rational(1, 2)
count = 0


def circle(expr):
    """Rational equality in R[L,r,u,c,s,...]/(s^2+c^2-1)."""
    numerator, denominator = S.fraction(S.factor(expr))
    remainder = S.rem(S.Poly(numerator, s), S.Poly(s*s+c*c-1, s))
    return S.factor(remainder.as_expr()/denominator)


def check(label, expr, on_circle=False):
    global count
    value = circle(expr) if on_circle else S.factor(expr)
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    count += 1
    print(f"PASS {label}")


C0 = 1/(2*L*(1-r))
A0 = (1+L*(1-r))/(2*L*r)
a0 = L/(2*(1-L*r))
B0 = (1+L*(1-r))/(2*L*r*(1-r))
gamma0 = (1+L*(1-r))/(2*L)
check("leading source B", C0+A0+h-B0)
check("leading source ratio", (C0+h)/B0-r)
check("leading source gamma", (C0+h)*A0/B0-gamma0)
check("leading gamma alternative", gamma0-h-1/(4*a0))

M0 = (C0+c*c*h)*(P-Y)**2+(A0+s*s*h)*Y**2+a0*(s*Y-T)**2-s*P*T
M0_sos = a0*(T-s*Y-(1-L*r)*s*P/L)**2+B0*(Y-r*P)**2+gamma0*(c*P)**2
check("three-square compact Bob completion", M0-M0_sos, True)

features6 = S.Matrix([T, s*Y, s*P, Y, P, c*P])
G6 = S.Matrix([
    [a0, -a0, -h, 0, 0, 0],
    [-a0, a0, h, 0, 0, 0],
    [-h, h, 1/(4*a0), 0, 0, 0],
    [0, 0, 0, B0, -B0*r, 0],
    [0, 0, 0, -B0*r, B0*r*r, 0],
    [0, 0, 0, 0, 0, gamma0],
])
check("evaluated six-feature Gram", (features6.T*G6*features6)[0]-M0_sos)

K = (X*X+c*c*(P-Y)**2+s*s*Y*Y-2*s*P*T
     +2*u*c*X*(s*P-s*Y+T)+u*u*s*s*X*X+u*u*(s*Y-T)**2)
K_original = ((c*(P-Y)+u*s*X-s*Z)**2
              +(1+u*u)*(s*Y-c*Z)**2+(X+u*Z)**2-(1+u*u)*Z*Z)
check("complete original core", K_original-K.subs(T, c*Z), True)

def selected(Cv, Av, av):
    return Cv*(P-Y)**2+Av*(1+u*u)*Y*Y+av*(1+u*u)*(s*Y-T)**2+K/2

S0 = selected(C0, A0, a0)
f = S.symbols("f")
event_excesses = (C-f*C0)*(P-Y)**2+(A-f*A0)*(1+u*u)*Y*Y+(a-f*a0)*(1+u*u)*(s*Y-T)**2
check("actual-source common-factor reassembly", selected(C,A,a)-(f*S0+(1-f)*K/2+event_excesses))
raw = S.Matrix([P, X, Y, T])
Graw = S.hessian(S0, (P, X, Y, T))/2
restricted_sos = (C0*(P-Y)**2+A0*(1+u*u)*Y*Y+a0*(1+u*u)*s*s*Y*Y
                  +(c*(P-Y)+u*s*X)**2/2+X*X/2+s*s*(1+u*u)*Y*Y/2)
check("strict full upper state block", S0.subs(T, 0)-restricted_sos)
det = ((1+u*u)**2*(L*c*c+u*u*s*s)*(1+L*(1-r))**2
       /(16*L*L*r*(1-r)*(1-L*r)))
check("full four-state determinant", Graw.det(method="domain-ge")-det, True)

t, q, e = u*u, c*c, 1-r
N2 = 1+r*r+e*q*(1+r*(1-q))+e*(2-q*r)*t+e*e*(1-q)*t*t
N1 = 2*e*q+(2-r*q*(3-q))*t+e*(1-q)*t*t
N0 = 1+(2-q)*t
trace_formula = 2*(L*L*N2+L*N1+N0)/((1+t)*(L*q+t*(1-q))*(1+L*e))
trace_adjugate = sum(Graw.minor_submatrix(i, i).det(method="domain-ge") for i in range(4))
check("opening-cleared trace inverse", trace_adjugate-det*trace_formula, True)

# Complete event-coordinate Gram, with no test/source residue choice.
def unit(i, j):
    M = S.zeros(4)
    M[i-1, j-1] = 1
    return M

eta = [u*u, c*c, s*s, s, u*c*s, u*c, u*u*s*s]
As = [
    S.diag(0, Ap, ap+h, 0),
    unit(1, 1)/2,
    -(unit(1, 2)+unit(2, 1)+unit(2, 2))/2,
    (unit(1, 3)+unit(3, 1)+unit(2, 3)+unit(3, 2))/2,
    (unit(1, 4)+unit(4, 1)+unit(2, 4)+unit(4, 2))/2,
    -(unit(3, 4)+unit(4, 3))/2,
    unit(4, 4)/2,
]
Gconst = S.diag(Cp, Ap, ap, h)
Gevent = Gconst + sum((f*M for f, M in zip(eta, As)), S.zeros(4))
events = S.Matrix([P-Y, Y, s*Y-T, X])
check("all seven evaluated event-coordinate matrices",
      (events.T*Gevent*events)[0]-selected(Cp, Ap, ap))

m, D, beta = S.symbols("m D beta")
rm = s*W-c*X+u*s*Y-u*c*Z
rD = c*(u*W-Y)+u*s*X-s*Z
rbeta = X+u*Z
full_original = (m*rm**2+C*(u*W-Y)**2+A*(1+u*u)*Y*Y
                 +(D+h)*rD**2+(a+h)*(1+u*u)*(s*Y-c*Z)**2
                 +(beta+h)*rbeta**2-h*(1+u*u)*Z*Z)
full_restored = selected(C,A,a).subs({P:u*W,T:c*Z}, simultaneous=True)+m*rm**2+D*rD**2+beta*rbeta**2
check("complete C1 parent and both normal/event restorations", full_original-full_restored, True)

# The fixed polarization identity is checked with symbolic test increments.
delta, cap = S.symbols("delta Bstar")
b = S.symbols("b1:8")
Gcenter = S.Matrix(4, 4, lambda i, j: S.Symbol(f"g{i}{j}") if i <= j else S.Symbol(f"g{j}{i}"))
target = Gcenter + sum((bj*Mj for bj, Mj in zip(b, As)), S.zeros(4))
polarized = Gcenter-7*delta*cap*S.eye(4)
for bj, Mj in zip(b, As):
    polarized += ((cap*S.eye(4)+Mj)*(delta+bj)+(cap*S.eye(4)-Mj)*(delta-bj))/2
check("fixed-radical block identity", sum((target-polarized).applyfunc(S.expand)))
assert all(S.expand(entry) == 0 for entry in target-polarized)

k = C+h
B = C+A+h
gamma = h+1/(4*a)
g = k*A/B
Delta = 8*C*A*a-4*C*a-1-2*(C+A+a)
defect = Delta/(8*a*B)
check("defect coefficient", defect-g+gamma)
check("leading three-event wall", Delta.subs({C:C0, A:A0, a:a0}))
M = (C+c*c*h)*(P-Y)**2+(A+s*s*h)*Y*Y+a*(s*Y-T)**2-s*P*T
M_defect = a*(T-s*Y-s*P/(2*a))**2+B*(Y-k*P/B)**2+g*c*c*P*P+defect*s*s*P*P
check("general exact below-wall defect", M-M_defect, True)
F = (a+h)/(4*a*a)
E = A*(k/B)**2+F
witness_selected = selected(C, A, a).subs({P:1, X:0, Y:k/B, T:s*k/B+s/(2*a)}, simultaneous=True)
check("complete selected-state witness", witness_selected-(defect+gamma*c*c+E*u*u-F*u*u*c*c), True)
wt = -defect/(2*(gamma+E+1))
check("explicit finite negative value", defect+(gamma+E)*wt-F*wt*wt-(defect*(gamma+E+2)/(2*(gamma+E+1))-F*wt*wt))

# A concrete symbolic-free rational coefficient sample also checks the six-Gram.
sample = {L:S.Rational(1), r:S.Rational(1,3)}
G6sample = G6.subs(sample)
assert G6sample.rank() == 3
assert all(S.det(G6sample.extract(I,I)) >= 0
           for size in range(1,7)
           for I in __import__("itertools").combinations(range(6),size))
count += 1
print("PASS rational six-feature Gram example, all principal minors, rank three")

def entries(M):
    return [[str(S.factor(M[i,j])) for j in range(M.cols)] for i in range(M.rows)]

data = {
    "status": "leading compact-Bob completion, actual C1 fixed-annulus family, and exact defect; no full-Q222 adequacy claim",
    "base_commit": "fb58c232e374dadf3cbb5f23e7e7413d429de884",
    "leading_source_domain": "L>0, 0<r<1, L*r<1",
    "testing_circle": "c>=0, c^2+s^2=1",
    "raw_variables": ["P", "X", "Y", "T"],
    "physical_substitution": {"P":"u*W", "T":"c*Z", "c":"1/sqrt(1+v^2)", "s":"v*c"},
    "leading_weights": {"C0":str(C0), "A0":str(A0), "a0":str(a0)},
    "six_functions": [str(f) for f in features6],
    "six_gram": entries(G6),
    "six_gram_rational_example": {"L":"1", "r":"1/3", "matrix":entries(G6sample), "rank":3},
    "event_coordinates": [str(f) for f in events],
    "constant_matrix": entries(Gconst),
    "testing_monomials": [str(f) for f in eta],
    "coefficient_matrices": [entries(M) for M in As],
    "fixed_annulus": {
        "source_domain":"0<ell<=L<=U, 0<r<1, L*r<1, C,A,a>0",
        "fixed_constant_domain":"ell,U,tau are fixed positive rationals, ell<=U, tau<=1",
        "actual_source_admission":"abs(C/C0-1)+abs(A/A0-1)+abs(a/a0-1)<=kstar/(2*(kstar+3))",
        "testing_domain":"abs(u)<=1, c>=0, c^2+s^2=1, u^2+c^2>=tau^2, 0<tau<=1",
        "kstar":"min(ell,1)*tau^2/(8*(7*U^2+5*U+3))",
        "N":"4+216/kstar",
        "hstar":"min(1,kstar/24)",
        "Bstar":"N+1",
        "n":"ceil(max(8/tau^2,56*Bstar/hstar))",
        "delta":"4/n",
        "centers":"(j/n,cos(k*pi/(2*n)),sin(k*pi/(2*n))), -n<=j,k<=n, retain u_j^2+c_k^2>=tau^2/2",
        "charts":"circle points within 1/n of center in each of u,c,s",
        "main_functions":"e_i and e_i*sqrt(delta +/- (eta_j-eta_j(center))), i=1..4,j=1..7",
        "retained_functions":["sqrt(1+u^2)*Y", "sqrt(1+u^2)*(s*Y-T)"],
        "main_source_gram":"diag(Gprime(center)-7*delta*Bstar*I4,(Bstar*I4+A_j)/2,(Bstar*I4-A_j)/2 for j=1..7)",
        "source_clipping":"Cprime=min(C,N), Aprime=min(A,N), aprime=min(a,N)",
        "excess_source_gram":"add C-Cprime to e_1 diagonal; append A-Aprime and a-aprime",
        "selected_functions_per_chart_max":62,
        "full_C1_extra_functions":[str(rm),str(rD),str(rbeta)],
        "full_C1_extra_source_coefficients":["m","D","beta"],
        "full_C1_functions_per_chart_max":65,
        "charts_max":"(2*n+1)^2",
    },
    "below_wall": {
        "Delta":str(Delta),
        "defect":str(defect),
        "witness":"P=1,X=0,Y=k/B,T=s*(k/B+1/(2*a)),c^2=u^2=-defect/(2*(gamma+E+1))",
        "witness_value":"defect*(gamma+E+2)/(2*(gamma+E+1))-F*t^2, F=(a+1/2)/(4*a^2), t=c^2=u^2",
    },
    "exact_checks_passed":count,
}
dest = Path(__file__).with_name("compact_bob_matrices.json")
dest.write_text(json.dumps(data, indent=2)+"\n")
print(f"All {count} exact checks passed; wrote {dest.name}.")
