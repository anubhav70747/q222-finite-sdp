"""Portable exact structural checks; Python 3 and SymPy only.

This checks printed algebra, coefficient-majorant arithmetic, features,
and rational margins. Source inequalities are proved in the note and
independent source review, not inferred from this output.
"""
import json
from pathlib import Path
import sympy as s

checks = {}

def zero(name, expr):
    checks[name] = s.expand(expr) == 0
    assert checks[name], name

u, w, x, sig, q, L, H0 = s.symbols('u w x sig q L H0')
W, X, Y, Z, b, h, e = s.symbols('W X Y Z b h e')
m, C, D, beta, a, ell, p, z = s.symbols('m C D beta a ell p z')
kap, H, dG, nu, R0, Rx, Bx, Ay, chi = s.symbols(
    'kap H dG nu R0 Rx Bx Ay chi')
Eyy, ERR, EyR = s.symbols('Eyy ERR EyR')
N, JC, JL, T, V, A = s.symbols('N JC JL T V A')
Gamma = T-L*V/2
R = JC-JL
Ny = (1-b)*JC+b*JL
Kc = 1+L**2

zero('source_baseline', 2*H0*Gamma**2+L*Gamma*V+Kc*V**2
     -(T**2+V**2+(2*H0-1)*Gamma**2+3*L**2*V**2/4))

Ac = q*(u**2-x)*W
Vc = u*w*(u**2-x)*W
Gc = sig*Z-u*w*W
Rc = w*Z/q-u*q*W
zero('zero_correction', Vc**2-sig*u**2*Ac**2
     -u*sig*(u**2-x)*Ac*Rc+(u**2-x)*Vc*Gc)

Eplus = (4*beta*(1+w**2)*N**2+2*(2*D+1)*s.Symbol('uXe')**2
         +4*b*p*z*u**2*Gamma**2
         +4*a*(u**2+w**2+u**2*w**2)*JC**2
         +4*ell*u**2*JL**2+4*C*z*w**2*JC**2)
raw = (Eplus+(2*(1+w**2)+4*m*w**2)*N**2
       +(4*a+2*z*(2*C+1))*JC**2+4*ell*JL**2
       +Eyy*Ny**2+ERR*R**2+EyR*Ny*R+dG*Gamma**2
       -4*nu*u*N*Gamma+q*w*(R0-8*p*u**2)*JC*Gamma
       -4*u*q*w*(2*m+2*D+1)*N*Ny+u*q*w*Rx*N*R
       +q*w*Bx*N*A+u*Ay*Ny*A+chi*u*A*R
       +4*kap*H*A**2+L*u**2*Gamma**2+L*Gamma*V+sig*L*u*A*R)
repair = Kc*(V**2-sig*u**2*A**2-u*sig*(u**2-x)*A*R
             +(u**2-x)*V*Gamma)
correction = (Eyy*Ny**2+ERR*R**2+EyR*Ny*R
              +(dG-2*H0)*Gamma**2-4*nu*u*N*Gamma
              +q*w*(R0-8*p*u**2)*JC*Gamma
              -4*u*q*w*(2*m+2*D+1)*N*Ny+u*q*w*Rx*N*R
              +q*w*Bx*N*A+u*Ay*Ny*A
              +(chi+sig*L)*u*A*R+(4*kap*H-2*H0)*A**2
              +L*u**2*Gamma**2-Kc*sig*u**2*A**2
              -Kc*sig*u*(u**2-x)*A*R+Kc*(u**2-x)*V*Gamma)
retained = (4*beta*N**2+(4*beta+2+4*m)*w**2*N**2
            +2*(2*D+1)*s.Symbol('uXe')**2
            +4*b*p*z*u**2*Gamma**2
            +4*a*(u**2+w**2+u**2*w**2)*JC**2
            +4*ell*u**2*JL**2+4*C*z*w**2*JC**2
            +(4*a+2*z*(2*C+1)-1)*JC**2+(4*ell-1)*JL**2
            +(2*H0-1)*Gamma**2+3*L**2*V**2/4
            +2*(H0-1)*A**2)
diagonal = 2*N**2+JC**2+JL**2+T**2+V**2+2*A**2
zero('complete_correction_from_SP6', raw+repair-retained-diagonal-correction)

eps = s.Symbol('eps', nonnegative=True)
n, Cstar = JC+JL, T+4*V
parts = [(4*u**2+4*w**2)*n**2,
         (64*eps+16*u**2+4*w**2)*n**2,
         (32*eps+12*u**2+8*w**2)*n**2,
         200*eps*Cstar**2, 12*u*N*Cstar,
         (512*eps*w+16*u**2*w)*JC*Cstar,
         16*u*w*N*n,32*u*w*N*n,48*eps*w*N*A,
         32*u*n*A,65*u*A*n,400*eps*A**2,8*u**2*Cstar**2,
         65*eps*u**2*A**2,65*eps*(u**3+u)*A*n,
         65*(u**2+eps)*V*Cstar]
roots = [N,JC,JL,T,V,A]
one = {v:1 for v in [u,w,eps]+roots}
norms = [int(part.subs(one)) for part in parts]
Qmajor = s.expand(sum(parts))
constant = s.expand(Qmajor.subs({u:0,w:0}))
assert sum(norms) == 10189
assert s.expand(constant.subs({r:1 for r in roots})) == 6109*eps
support = set(s.Poly(Qmajor,u,w).monoms())
assert support == {(0,0),(0,1),(0,2),(1,0),(1,1),(2,0),(2,1),(3,0)}
checks.update(majorant_norm=sum(norms),constant_norm=6109,
              majorant_parts=norms,testing_support=sorted(support))

epsv, deltav = s.Rational(1,2**160),s.Rational(1,2**100)
c0, Astar = s.Rational(1,256),s.Integer(2**30)
gam = s.Rational(1,7*2**38)
assert s.Rational(11,16)-2**20*epsv > c0
assert s.Rational(1,32)-2**19*epsv > c0
assert deltav < gam and 7*Astar*gam == c0
assert 2-200*epsv > 0
checks['rational_margins'] = True
q0, aj, bj = s.symbols('q0 aj bj')
zero('scalar_channel_identity',
     (q0-c0)+6*Astar*gam+((Astar+aj)*(gam+bj)
       +(Astar-aj)*(gam-bj))/2-(q0+aj*bj))

F0 = [X,Y,Z,W,u*Z,w*Z,u*W,u*w*W,u**2*W,u**3*w*W]
extra = [w*X,u*w*Z,u*X,u**2*w*W,u*Y,w*Y,u*w*Y]
allowed = {s.Poly(v,u,w,W,X,Y,Z).monoms()[0] for v in F0+extra}
Nxc, Jcc, Jlc = X+u*h*Z,Y/q-b*u*q*W,Y/q-w*Z/q+(1-b)*u*q*W
Tc = Gc+L*Vc/2
root_expr = [Nxc,Jcc,Jlc,Tc,Vc,Ac]
ret_expr = [Nxc,w*Nxc,u*(X+e*w*W),u*Gc,u*Jcc,w*Jcc,
            u*w*Jcc,u*Jlc,Jcc,Jlc,Gc,Vc,Ac]
base_allowed = {s.Poly(v,u,w,W,X,Y,Z).monoms()[0] for v in F0}
for expr in root_expr:
    assert set(s.Poly(s.expand(expr),u,w,W,X,Y,Z).monoms()) <= base_allowed
for expr in ret_expr:
    assert set(s.Poly(s.expand(expr),u,w,W,X,Y,Z).monoms()) <= allowed
checks['feature_counts'] = dict(base=10,main=150,retained_additional=7,total=157)
checks['fixed_feature_spans'] = True

# Nonempty source path, uniformly for16<=B<=18 and0<t<=epsilon.
t, B = s.symbols('t B')
zp=bp=1-t**2;rp=bp*zp;hp=t-4*t**2+B*t**3;thp=1-hp
Jp=1+t*hp-t*(1-rp)
Pp=s.expand(t+1+t+t*zp*(1+bp*t+2*t*hp)
            +t*bp*t**2*hp*zp**2+t**2*hp*thp*zp)
pNp=s.expand(zp*(1+t*hp)*(1-hp)+(1-rp)*(1+t*hp*zp))
DNp=s.expand((1+t)*(1+t*zp)*(1+t*rp)-Jp*Pp)
bNp=s.expand((1+t)*(1+t*rp)*(1+t*hp*zp)-thp*Pp)
LNp=s.expand(4*DNp*hp+2*bNp*thp*Jp+2*hp**2*Jp*Pp
              -2*bp*pNp*t**2*zp*Jp)
path_norms=[]
for expr,power,leading,expected in [
    (DNp,2,1,27640),(bNp,3,B-17,24459),(LNp,3,2*B-38,248856)]:
    assert s.Poly(expr,t).coeff_monomial(t**power)==leading
    norm=sum(abs(c)*18**j for (i,j),c in s.Poly(expr,t,B).terms() if i>power)
    assert norm==expected
    path_norms.append(int(norm))
assert 27640*epsv<1 and 24459*epsv<1 and 248856*epsv<2
checks['nonempty_source_path_remainder_norms']=path_norms

# Actual local data, not only a feature count. All expressions use SymPy
# syntax; sqrt is the positive radical in the stated labelled field.
gs = h-(1-b)
k1 = 2*(1+2*m)*gs**2+8*D*h*b+8*b*p*sig+4*beta*(1-h)**2
coefficient_definitions = {
    'Eyy':2*z*((1+2*m)*u**2+(1+2*D)*w**2),
    'ERR':z*(8*b*D+8*b*p*sig+4*b*sig*(2*C+2*D+1)*(b*z-1)
             +k1*u**2+2*b**2*(2*D+1)*w**2),
    'EyR':4*z*(u**2*(2*D*h+2*p*sig+(1+2*m)*gs)
             +b*(2*D+1)*w**2+sig*(2*C+1)*(b*z-1)
             +2*D*(1+sig*(b*z-1))+2*p*sig)}
quad = diagonal+correction.subs({Eyy:coefficient_definitions['Eyy'],
                               ERR:coefficient_definitions['ERR'],
                               EyR:coefficient_definitions['EyR']})
G = s.hessian(quad,roots)/2
monomials = [w,w**2,u,u*w,u**2,u**2*w,u**3]
base_monomials = [s.Poly(f,u,w,W,X,Y,Z).monoms()[0] for f in F0]
all_monomials = [s.Poly(f,u,w,W,X,Y,Z).monoms()[0] for f in F0+extra]
def vector(expr, mons):
    poly = s.Poly(s.expand(expr),u,w,W,X,Y,Z)
    return [s.sstr(poly.coeff_monomial(mon)) for mon in mons]

root_rows = [vector(expr,base_monomials) for expr in root_expr]
directions=[]
for i in range(6):
    d=[0]*6;d[i]=1
    directions.append((f'diagonal_{i+1}',d,G[i,i]-s.Rational(5,16)))
for i in range(6):
    for j in range(i+1,6):
        for sign in [1,-1]:
            d=[0]*6;d[i]=1;d[j]=sign
            directions.append((f'pair_{i+1}_{j+1}_{"plus" if sign==1 else "minus"}',
                               d,(s.Rational(1,16)+sign*G[i,j])/2))
channel_data=[]
for name,d,expr in directions:
    poly=s.Poly(s.expand(expr),u,w)
    channel_data.append({'id':name,'fixed_root_direction':d,
        'constant_source_coefficient':s.sstr(poly.coeff_monomial(1)),
        'nonconstant_source_coefficients':[
             s.sstr(poly.coeff_monomial(mon)) for mon in monomials]})
assert len(channel_data)==36

features = [s.sstr(f) for f in F0]
for mon in monomials:
    for f in F0:
        features += [s.sstr(f*s.sqrt(gam+mon)),s.sstr(f*s.sqrt(gam-mon))]
features += [s.sstr(f) for f in extra]
assert len(features)==157 and len(set(features))==157

retained_data=[]
for coef,expr in [
    (4*beta,Nxc),(4*beta+2+4*m,w*Nxc),
    (2*(2*D+1),u*(X+e*w*W)),(4*b*p*z,u*Gc),
    (4*a,u*Jcc),(4*a,w*Jcc),(4*a,u*w*Jcc),
    (4*ell,u*Jlc),(4*C*z,w*Jcc),
    (4*a+2*z*(2*C+1)-1,Jcc),(4*ell-1,Jlc),
    (2*H0-1,Gc),(3*L**2/4,Vc),(2*(H0-1),Ac)]:
    retained_data.append({'nonnegative_source_coefficient':s.sstr(coef),
         'root_coefficients_on_17_polynomial_features':vector(expr,all_monomials)})

data = {
 'scope':'Complete local lambda<=0 C1 parent on the stated source/testing box; not a global pencil.',
 'symbol_conventions':{'sig':'sigma','kap':'kappa','ell':'opening-cleared alpha*z',
   'q':'sqrt(z), source only','N':'Nx','JC':'J_C','JL':'J_ell',
   'T':'Gamma+(L/2)*V','A':'A_theta'},
 'epsilon':str(epsv),'delta':str(deltav),'gamma':str(gam),
 'c0':str(c0),'Astar':str(Astar),
 'source_conditions':['0<x,sigma<=epsilon','0<z<=1','0<b<1',
   'r=b*z<1','0<h<1','E=1-x*h/sigma>0','native J,M>0',
   'D,beta>=0','L=lambda/sigma<=0'],
 'testing_conditions':['0<=u,w<=delta','W,X,Y,Z unrestricted real'],
 'source_definitions':{
   'r':'b*z','Q':'1-b','theta':'1-h','E':'1-x*h/sig',
   'K':'1+x','Rs':'1+x*r','M':'Rs-x*z*theta',
   'J':'1+x*h-sig*(1-r)',
   'P':'sig+K+sig*z*(1+b*x+2*x*h)+sig*b*x**2*h*z**2+x**2*h*theta*z',
   'p':'(z*(1+x*h)*E+(1-r)*(1+x*h*z))/(2*P)',
   'l':'(1+sig*r)*(e*Rs*(1+sig*z)+M*E)/(2*P)',
   'kap':'K*(1+sig*r)*(1+sig*z)/(2*M*P)',
   'm':'kap*x**2*z*H',
   'C':'K*z*(1+sig*r)*E/(2*(1-r)*P)',
   'D':'(K*(1+sig*z)*(1+sig*r)-J*P)/(2*J*P)',
   'beta':'(K*(1+sig*r)*(1+x*h*z)-theta*P)/(2*theta*P)',
   'q':'sqrt(z)','H':'h/sig','H0':'1+z*(1-b)','e':'x*h/sig',
   'dG':'4*D/sig**2','L':'lambda/sig','a':'p/b','ell':'l/(1-b)',
   'event_A':'a/z','event_alpha':'ell/z','k':'sqrt(x)',
   'nu':'(2*D-2*beta*(1-h)+h)/sig',
   'lambda':'(8*D*h+4*beta*(1-h)**2+2*h**2)/sig-4*b*p*sig*z',
   'R0':'-4*(2*(C+D)*(b*z-1)+2*p+b*z-1)',
   'Rx':s.sstr(-4*(2*D*b+(1+2*m)*gs-2*beta*(1-h))),
   'Bx':'8*kap*e',
   'Ay':s.sstr(8*(p-s.Symbol('l'))+4*z*((1+2*m)*gs+2*D*h+2*p*sig)),
   'chi':s.sstr(8*(p*b+s.Symbol('l')*(1-b))-4*h+4*e*h*(2*D+1)-8*D*h
        +4*z*(4*D*b*h+4*p*b*sig+(1+2*m)*gs**2+2*beta*(1-h)**2))},
 'shared_native_source_definition_file':'fixed_box_c1.tex, section1',
 'root_order':['Nx','JC','Jl','T','V','Atheta'],
 'root_coefficients_on_10_polynomial_features':root_rows,
 'polynomial_features_10':[s.sstr(f) for f in F0],
 'retained_additional_polynomial_features_7':[s.sstr(f) for f in extra],
 'fixed_testing_features_157':[{'id':f'f{i+1:03d}','expression':expr}
                             for i,expr in enumerate(features)],
 'testing_channel_monomials_7':[s.sstr(mon) for mon in monomials],
 'gram_entries_6_by_6':[[s.sstr(G[i,j]) for j in range(6)] for i in range(6)],
 'main_channels_36':channel_data,
 'retained_squares_14':retained_data,
 'source_bounds':{'H':2,'abs_L':8,'abs_R0':'512*epsilon',
    'gram_constant_entry_error':'2**20*epsilon',
    'gram_nonconstant_coefficient_bound':'2**20'},
 'channel_factorization':{
   'identity':'q=(q0-c0)+sum_j((Astar+a_j)*(gamma+b_j)+(Astar-a_j)*(gamma-b_j))/2',
   'base_multiplier_source_weight':'q0-c0',
   'plus_multiplier_source_weight':'(Astar+a_j)/2',
   'minus_multiplier_source_weight':'(Astar-a_j)/2',
   'direction_row':'fixed_root_direction * root_coefficients_on_10_polynomial_features'},
 'zero_faces':'beta=0 literal; D=0 cannot intersect this branch since dG>=2-200*epsilon>0',
 'nonempty_source_path':{'x':'t','sig':'t','z':'1-t**2','b':'1-t**2',
   'h':'t-4*t**2+B*t**3','parameter_bounds':'0<t<=epsilon,16<=B<=18',
   'strict_source':'B=18 gives D>0,beta>0,L<0',
   'beta_zero':'An algebraic B in(16,18) gives beta=0,D>0,L<0',
   'remainder_norms':path_norms},
 'signed_testing_gauge':{'u':'abs(us)','w':'abs(ws)','W':'Ws',
       'X':'sign_plus(ws)*Xs','Y':'sign_plus(us)*Ys',
       'Z':'sign_plus(us)*sign_plus(ws)*Zs','sign_plus(0)':1},
 'complex_states':'Apply the real certificate to real and imaginary parts separately.'}
outdir=Path(__file__).resolve().parent
(outdir/'fixed_box_c1_nonpositive_data.json').write_text(
    json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8')
(outdir/'fixed_box_c1_nonpositive_checks.json').write_text(
    json.dumps(checks,indent=2,sort_keys=True)+'\n',encoding='utf-8')
checks['data_file']='fixed_box_c1_nonpositive_data.json'
print(json.dumps(checks,indent=2))
