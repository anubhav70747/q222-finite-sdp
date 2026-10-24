"""Part A: source checks for R6.
1) native source (c2-native-source) is stationary at native contact for parent (c2-master-parent)
2) c2-ba-source == native source under the weak-II map
3) R6 substitution a=1/s, B=s^2 b, e=s(1+hsb-g): residues (RB3),(RB4)
"""
import sympy as sp

u, v, W, X, Y, Z = sp.symbols('u v W X Y Z', real=True)
n_, m_, C_, A_, D_, al_ = sp.symbols('n m C A D alpha', real=True)

def parent(n, m, C, A, D, al, u=u, v=v, W=W, X=X, Y=Y, Z=Z):
    Du, Dv = 1 + u**2, 1 + v**2
    M = v*W - X + u*v*Y - u*Z
    J = u*W + u*v*X - Y - v*Z
    return (n*Du*Dv*X**2 + m*M**2 + C*Dv*(u*W - Y)**2 + A*Du*Dv*Y**2
            + (D + sp.Rational(1, 2))*J**2 + (al + sp.Rational(1, 2))*Du*(v*Y - Z)**2
            + sp.Rational(1, 2)*Dv*(X + u*Z)**2 - sp.Rational(1, 2)*Du*Dv*Z**2)

# --- 1) native source
c, an, bn, tn, r = sp.symbols('c a_n b_n t_n r', positive=True)
z = c**2
Ac = 1 + z*r
hn = bn - an**2*tn
gn = (1 + z)*bn - (1 - r + an**2*z*tn)
jn = 1 + z*r - tn
Delta = ((1 + z)*z*bn**2 + 2*z*(1 + z)*bn*r + (1 + z)*(1 + an**2*z)*bn
         + z**2*(1 - an**2)*r**2 + 2*z*(1 + an**2*z)*r + 1 + an**2*z)
Nn = (1 + z)*(bn + r)*(bn + an**2*(1 - tn)) - an**2*(1 + an**2*z)*(1 - r)*tn
nat = dict(
    n=((bn + r)*Ac*(Ac + an**2*z*tn) - r*Delta)/(2*r*Delta),
    m=(1 + an**2*z)*(1 - r)*Ac*(Ac + an**2*z*tn)/(2*gn*Delta),
    C=Ac*Nn/(2*an**2*tn*Delta),
    A=an**2*Ac*(bn + r)*(Ac + an**2*z*tn)/(2*hn*Delta),
    D=(1 + z)*(bn + r)*(bn - an**2*z*r)*(Ac + an**2*z*tn)/(2*an**2*jn*Delta),
    al=z**2*(bn + r)*Nn/(2*(1 - z*hn)*Delta),
)
contact = {u: c, v: an*c, W: bn/an, X: -c*r, Y: c*(bn/an - an*tn), Z: 1}
Pn = parent(nat['n'], nat['m'], nat['C'], nat['A'], nat['D'], nat['al'])
vals = [Pn] + [sp.diff(Pn, w) for w in (u, v, W, X, Y, Z)]
res = [sp.simplify(sp.together(e.subs(contact))) for e in vals]
print('native stationarity (value, d/du, d/dv, d/dW, d/dX, d/dY, d/dZ):', res)

# --- 2) c2-ba-source vs native under weak-II map
a, B, t, h = sp.symbols('a B t h', positive=True)
l = 1 - h
e = 1 - r
x = a*B
g = 1 - a*e + h*a*B
j = 1 + a*B*r - l/t
ca = 1 + a*B*r
S = 1 + B*(l + a*r)
Dl = (a + (1 + a*B) + (1 + a*B)*B + 2*a*B*(1 + a*B)*r + a**3*B**2*r**2 + 2*a**2*B*r
      + t*B*(1 + a + a*B + a**2*B*r*(2 - r)))
Mb = (1 + a*B)*(1 + a*r)*(h + t) - a*(1 + t*B)*e*l
ba = dict(
    n=((1 + a*r)*ca*S - r*Dl)/(2*r*Dl),
    m=a**2*(1 + t*B)*e*ca*S/(2*g*Dl),
    C=ca*Mb/(2*l*Dl),
    A=(t*ca*(1 + a*r)*S/(2*Dl))/h,
    D=(1 + a*B)*(1 + a*r)*(1 - t*a*B*r)*S/(2*t*j*Dl),
    al=B**2*(1 + a*r)*Mb/(2*(1 - h*B)*Dl),
)
# map: z=c^2=aB, an^2=t/a, bn=1/a, tn=l/t
mp = {c: sp.sqrt(a*B), an: sp.sqrt(t/a), bn: 1/a, tn: l/t}
for k in ba:
    d = sp.simplify(sp.together(nat[k].subs(mp) - ba[k]))
    print('ba-vs-native', k, d)

# --- 3) R6 substitution
s, b, gg = sp.symbols('s b g', positive=True)
sub = {a: 1/s, B: s**2*b}
# e = s(1+hsb-g)  <=> r = 1 - s(1+h s b - g)
sub_r = {r: 1 - s*(1 + h*s*b - gg)}
def R6(expr):
    return sp.factor(sp.together(expr.subs(sub).subs(sub_r)))
# check g identity
print('g check:', sp.simplify(R6(g) - gg))
quant = dict(
    n_over_s2=R6(ba['n'])/s**2,
    mg=R6(ba['m'])*gg,
    C=R6(ba['C']), A=R6(ba['A']), D=R6(ba['D']),
    al_over_s3b2=R6(ba['al'])/(s**3*b**2),
)
for k, q in quant.items():
    q = sp.factor(sp.together(q))
    num, den = sp.fraction(q)
    den0 = sp.factor(den.subs({s: 0, gg: 0}))
    den00 = sp.factor(den.subs({s: 0, gg: 0, b: 0}))
    resid = sp.factor(sp.together(q.subs({gg: 0})).subs(s, 0)) if den0 != 0 else 'den vanishes'
    print(f'{k}: den at s=g=0 -> {den0};  at s=g=b=0 -> {den00};  residue -> {resid}')
