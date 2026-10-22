"""Exact endpoint identities and rational constants for bounded-unit-L.
Requires SymPy; no solver or floating-point operation is used.
"""
from pathlib import Path
import json
import sympy as s

p,v,w,mu=s.symbols('p v w mu')
W,X,Y,Z,B,R,H=s.symbols('W X Y Z B R H')
Dp,Dv,Dw=1+p*p,1+v*v,1+w*w
checks={}
def ck(name,expr):
    if s.factor(expr)!=0:
        raise AssertionError(name+': '+str(s.factor(expr)))
    checks[name]=True

Kp=B*B-2*p*v*B*W-2*v*B*Y+Dp*R*R+2*Dp*R*W+(W-p*Y)**2+Dp*v*v*Y*Y
Kpw=H*H-2*p*H*W+2*p*w*H*X-2*H*Y+w*w*W*W+2*w*W*X-2*p*w*w*W*Y+(Dp+p*p*w*w)*(X*X+Y*Y)-2*p*w*X*Y
Qp=mu*Dv*(W-p*Y)**2+mu*Dp*Dv*Y*Y+Kp/2
Qpw=mu*Dw*(W-p*Y)**2+mu*Dp*Dw*Y*Y+Kpw/2
Gp=s.Matrix([
 [mu*Dv+s.Rational(1,2),-p*(mu*Dv+s.Rational(1,2)),-p*v/2,Dp/2],
 [-p*(mu*Dv+s.Rational(1,2)),mu*Dv*(1+2*p*p)+(p*p+Dp*v*v)/2,-v/2,0],
 [-p*v/2,-v/2,s.Rational(1,2),0],
 [Dp/2,0,0,Dp/2]])
Gpw=s.Matrix([
 [mu*Dw+w*w/2,-p*(mu*Dw+w*w/2),w/2,-p/2],
 [-p*(mu*Dw+w*w/2),mu*Dw*(1+2*p*p)+(1+p*p*Dw)/2,-p*w/2,-s.Rational(1,2)],
 [w/2,-p*w/2,(1+p*p*Dw)/2,p*w/2],
 [-p/2,-s.Rational(1,2),p*w/2,s.Rational(1,2)]])
ap=s.Matrix([W,Y,B,R]); aw=s.Matrix([W,Y,X,H])
ck('CE4 endpoint full selected matrix',(ap.T*Gp*ap)[0]-Qp)
ck('CE7 endpoint full selected matrix',(aw.T*Gpw*aw)[0]-Qpw)
ck('CE4 zero-p reserve squares',Qp.subs(p,0)-(mu*Dv*(W*W+Y*Y)+((W+R)**2+(B-v*Y)**2)/2))
ck('CE7 zero-p reserve squares',Qpw.subs(p,0)-(mu*Dw*(W*W+Y*Y)+((X+w*W)**2+(H-Y)**2)/2))
# Full raw quarter-parents in ordinary and reciprocal Bob coordinates.
C,A,alpha,m,D,beta=s.symbols('C A alpha m D beta')
corep=(W+v*X-p*Y-p*v*Z)**2+Dp*(v*Y-Z)**2+Dv*(p*X+Z)**2-Dp*Dv*Z*Z
fullp=(m*(p*v*W-p*X+v*Y-Z)**2+C*Dv*(W-p*Y)**2+A*Dp*Dv*Y*Y
       +(D+s.Rational(1,2))*(W+v*X-p*Y-p*v*Z)**2
       +(alpha+s.Rational(1,2))*Dp*(v*Y-Z)**2
       +(beta+s.Rational(1,2))*Dv*(p*X+Z)**2-Dp*Dv*Z*Z/2)
ck('CE4 full core restoration',Kp.subs({B:p*X+Z,R:v*X},simultaneous=True)-corep)
excessp=((C-mu)*Dv*(W-p*Y)**2+(A-mu)*Dp*Dv*Y*Y+alpha*Dp*(v*Y-Z)**2
         +m*(p*v*W-p*X+v*Y-Z)**2+D*(W+v*X-p*Y-p*v*Z)**2+beta*Dv*(p*X+Z)**2)
ck('CE4 all six original contributions',fullp-Qp.subs({B:p*X+Z,R:v*X},simultaneous=True)-excessp)
corepw=(w*W+X-p*w*Y-p*Z)**2+Dp*(Y-w*Z)**2+Dw*(p*X+Z)**2-Dp*Dw*Z*Z
ck('CE7 full core restoration',Kpw.subs(H,w*Z)-corepw)
fullpw=(m*(p*W-p*w*X+Y-w*Z)**2+C*Dw*(W-p*Y)**2+A*Dp*Dw*Y*Y
       +(D+s.Rational(1,2))*(w*W+X-p*w*Y-p*Z)**2
       +(alpha+s.Rational(1,2))*Dp*(Y-w*Z)**2
       +(beta+s.Rational(1,2))*Dw*(p*X+Z)**2-Dp*Dw*Z*Z/2)
excesspw=((C-mu)*Dw*(W-p*Y)**2+(A-mu)*Dp*Dw*Y*Y+alpha*Dp*(Y-w*Z)**2
         +m*(p*W-p*w*X+Y-w*Z)**2+D*(w*W+X-p*w*Y-p*Z)**2+beta*Dw*(p*X+Z)**2)
ck('CE7 all six original contributions',fullpw-Qpw.subs(H,w*Z)-excesspw)
# Exact termwise row-bound majorants for 0<=p<=1, mu<=1/4.
rowp=[p*(2*mu+1)+p*p/2,(4*mu+1)*p*p+p*(2*mu+s.Rational(1,2)),p/2,p*p]
rowpw=[p*(2*mu+1),(4*mu+1)*p*p+p*(2*mu+1),p+p*p,p]
def majorant(poly):
    result=0
    for (ip,im),co in s.Poly(s.expand(poly),p,mu).terms():
        assert co>=0 and ip>=1
        result+=co*p*s.Rational(1,4)**im
    return result
assert all(s.Poly(3*p-majorant(t),p).coeff_monomial(p)>=0 for t in rowp)
assert all(s.Poly(4*p-majorant(t),p).coeff_monomial(p)>=0 for t in rowpw)
checks['endpoint row-bound majorants']=True
assert s.Rational(1,3)-s.Rational(4,24)==s.Rational(1,6)
checks['endpoint positive reserve']=True
q,t=s.symbols('q t')
ck('generalized annulus denominator',(1+t)*(q+t*(1-q))-(q+t)-t*t*(1-q))
# Rational grid replacing trigonometric centers without changing counts.
tt=s.symbols('tt')
cc=(1-tt*tt)/(1+tt*tt); ss=2*tt/(1+tt*tt)
ck('rational compact circle grid',cc*cc+ss*ss-1)
ck('circle c derivative majorant',2*(1+tt*tt)**2-4*tt-(2*tt**4+4*(tt-s.Rational(1,2))**2+1))
ck('circle s derivative majorant',2*(1+tt*tt)**2-2*(1-tt*tt)-(2*tt**4+6*tt*tt))
# A completely evaluated rational instance U=1.
UU=s.Integer(1); MM=s.Integer(96); tau=s.Rational(1,2**25)
BB=UU**2*(4+2*MM**2+MM**4)+UU*(2+2*MM**2+MM**4)+1+2*MM**2
kk=tau**2/(8*BB); KK=3*MM**2
eps=min(s.Rational(1,2**78),kk/(50*(kk+KK)))
NN=KK+1+24*KK*KK/kk; hh=min(1,kk/24); bstar=NN+1
nn=s.ceiling(max(8*MM**2/tau**2,56*MM**2*bstar/hh))
J=(2*MM*nn+1)*(2*nn+1)
assert 25*eps<=kk/(2*(kk+KK))
assert 7*(4*MM**2/nn)*bstar<=hh/2
assert 4*MM**2/nn<=tau**2/2
checks['U1 all rational admission and grid inequalities']=True
def entries(M):return [[str(s.factor(a)) for a in row] for row in M.tolist()]
report={'scope':'bounded C1 source interval L in [ell,U], all tests; no global Q222 assertion',
        'checks':checks,'exact_checks':len(checks),'endpoint_matrices':{'p_v':entries(Gp),'p_w':entries(Gpw)},
        'endpoint_feature_count_per_Bob_chart':10,
        'U1_variant':'ell=U=1, optimized published-local669 option; tau=2^-25',
        'U1_rational_constants':{k:str(v) for k,v in {'U':UU,'M':MM,'tau':tau,'B':BB,'k':kk,'K':KK,'epsilon':eps,'N':NN,'h':hh,'Bstar':bstar,'n':nn,'chart_count_bound':J,'feature_count_bound_using_669_plus_annulus_plus_endpoint':669+65*J+20}.items()},
        'verification_limit':'The accompanying analytical proof establishes inequalities. Bounded-R inner admission is retained for general ell; local669 suffices for the optimized ell>=1 option.'}
Path(__file__).with_name('bounded_unit_L_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'exact_checks':len(checks),'U1_rational_constants':report['U1_rational_constants']},indent=2))
