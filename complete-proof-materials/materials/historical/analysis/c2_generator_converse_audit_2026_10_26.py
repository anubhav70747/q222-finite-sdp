"""Exact converse identities for positive ME C2 generators; no optimizer."""
import json
from pathlib import Path
import sympy as s

# Generic rational elimination in the converse, before imposing E=0.
a,tau,N,m,D=s.symbols('a tau N m D')
chi=1+(1-a*a)*tau
ri=(chi*N-a*m)/(N+m)
si=(chi*N-a*D)/(N+D)
E=((a+1)*chi*tau*N*N+(1-a*a)*tau*(1+tau)*(m+D)*N
   -(a+1)*(1+tau)*m*D)
den=(N+m)*(N+D)
gap=(1-ri)*(1-si)-(a+1)*tau*(a+ri+si-1)
assert s.cancel(gap+(a+1)*E/den)==0
assert s.cancel((chi-ri)*N-(a+ri)*m)==0
assert s.cancel((chi-si)*N-(a+si)*D)==0
assert s.cancel((chi-ri)*(1-si)-(a+1)*tau*(a*si+ri)+(a+1)*E/den)==0
assert s.cancel((chi-si)*(1-ri)-(a+1)*tau*(a*ri+si)+(a+1)*E/den)==0
assert s.cancel(1-ri-(a+1)*(m-(1-a)*tau*N)/(N+m))==0
assert s.cancel(1-si-(a+1)*(D-(1-a)*tau*N)/(N+D))==0
S=(m+D)/tau; k=a*a*tau/chi; B=k*S
P=(a+1)*N; d=2*(P-B); U=D/tau; V=m/tau; M=m+D
quartic=S*S*P*(M+d/2)-(a+1)**2*U*V*M*(B+S)
assert s.cancel(quartic-(a+1)*S*S*E/(chi*tau))==0

# The current C2_CONTACT_GEOMETRY.md dictionary, with d=1.
r,q=s.symbols('r s')
e=r+q; p=r*q; g=a+r+q-1
L=a*(e-p)+p; Q=a*g+L
Delta=a**3*p+a*a*(e*e-p)+a*p*(2*e-p)+p*p
B0=a**3*p+a*a*(e*e+p-e)+a*((3*e-2)*p-p*p)+p*p
nn=-a*B0/(2*(a+1)*Delta)
mm=(a+q)*(1-r)*(a*q+r)*L/(2*(a+1)*g*Delta)
CC=r*r*(a+q)*Q/(2*(1-r)*Delta)
AA=a*(a+r)*(a+q)*L/(2*(a+1)*Delta)
DD=(a+r)*(1-q)*(a*r+q)*L/(2*(a+1)*g*Delta)
al=q*q*(a+r)*Q/(2*(1-q)*Delta)
BB=AA+nn; UU=CC+s.Rational(1,2); VV=al+s.Rational(1,2)
SS=UU+VV; MM=mm+DD; PP=BB+s.Rational(1,2); NN=nn+s.Rational(1,2)
ta=s.cancel(MM/SS); ch=s.cancel(1+(1-a*a)*ta)
assert s.cancel(a*(1+2*nn)-2*AA)==0
assert s.cancel(ta-(1-r)*(1-q)/((a+1)*g))==0
assert s.cancel(BB/SS-a*a*(1-r)*(1-q)/((a+1)*L))==0
assert s.cancel(PP/UU-(1-r)*(a+q)/(a*r+q))==0
assert s.cancel(PP/VV-(1-q)*(a+r)/(a*q+r))==0
assert s.cancel((ch*NN-a*mm)/(NN+mm)-r)==0
assert s.cancel((ch*NN-a*DD)/(NN+DD)-q)==0
Gamma=a*a*(BB+SS)-BB
assert s.cancel(mm*UU-DD*VV)==0
assert s.cancel(mm*Gamma-BB*VV)==0
assert s.cancel(DD*Gamma-BB*UU)==0

# Original (unrescaled) native t and all its strict chart gaps.
tn=(1-r)/a; zn=q/(a*r); bn=a
hn=bn-a*a*tn
gn=(1+zn)*bn-(1-r+a*a*zn*tn)
jn=1+zn*r-tn
assert s.cancel(hn-a*r)==0
assert s.cancel(gn-g)==0
assert s.cancel(jn-g/a)==0
assert s.cancel(1-zn*hn-(1-q))==0
assert s.cancel(bn-a*a*zn*r-a*(1-q))==0

# a=1: the reduced quartic and minors yield the two Schur equalities.
At,mt,Dt=s.symbols('A m D')
Mt=mt+Dt; dt=4*At-2*Mt; nt=Mt-At
St=Mt*(At*At-mt*Dt)/(mt*Dt)
Ut=Dt*St/Mt; Vt=mt*St/Mt
assert s.cancel(Ut-dt/2-(At-mt)**2/mt)==0
assert s.cancel(Vt-dt/2-(At-Dt)**2/Dt)==0
assert s.cancel((At-mt)/(At+mt)-ri.subs({a:1,N:At,m:mt,D:Dt}))==0
assert s.cancel((At-Dt)/(At+Dt)-si.subs({a:1,N:At,m:mt,D:Dt}))==0

# The displayed nonphysical-branch rational source.
av=s.Rational(3,2); dv=s.Integer(1); nv=s.Rational(39,2)
mv=s.Rational(9,2); Cv=s.Rational(7,6); Av=s.Integer(30); Dv=s.Integer(1); alv=s.Integer(7)
Bv=Av+nv; Uv=Cv+dv/2; Vv=alv+dv/2; Sv=Uv+Vv; Mv=mv+Dv
Nv=nv+dv/2; Pv=Bv+dv/2; Tv=Mv+dv/2
tv=Mv/Sv; cv=1+(1-av*av)*tv
Gv=av*av*(Bv+Sv)-Bv
assert av*(dv+2*nv)==2*Av
assert mv*Uv==Dv*Vv and mv*Gv==Bv*Vv and Dv*Gv==Bv*Uv
assert Sv*Sv*Pv*(Mv+dv/2)==(av+1)**2*Uv*Vv*Mv*(Bv+Sv)==s.Rational(75625,3)
rv=(cv*Nv-av*mv)/(Nv+mv); sv=(cv*Nv-av*Dv)/(Nv+Dv)
assert rv==-s.Rational(1,14) and sv==s.Rational(1,6)
cosine=((Uv*Tv)**2*(Pv*Pv+Vv*Vv)-(Pv*Vv)**2*(Uv*Uv+Tv*Tv))/(2*Pv*Vv*Uv*Tv*(Pv*Vv+Uv*Tv))
assert cosine==-s.Rational(9,5)

# The regular finite r=0 endpoint in the response has no singular denominator.
endpoint=[a*(1-q)/(2*(a+1)*q),
 (a+q)/(2*(a+1)*(a+q-1)),0,
 a*(a+q)/(2*(a+1)*q),
 (1-q)/(2*(a+1)*(a+q-1)),(a+2*q-1)/(2*(1-q))]
assert all(s.cancel(expr.subs(r,0)-want)==0
           for expr,want in zip([nn,mm,CC,AA,DD,al],endpoint))

data={
 'main_generic_scalar_identity_and_quartic_reduction_exact':True,
 'all_native_ratio_minors_and_rational_inverse_identities_exact':True,
 'branch_inequalities_equivalent_to_0_less_r_s_less_1':True,
 'original_native_parameter_correction':'t_source=(1-r)/a, not 1-r; alternatively define rescaled a*t_source=1-r',
 'original_native_strict_gaps':['a*r','a+r+s-1','(a+r+s-1)/a','1-s','a*(1-s)'],
 'undefined_stationary_T_correction':'T=M+d/2 (named R in C2_CONTACT_GEOMETRY.md)',
 'a_equals_1_reduction_exact':True,
 'rational_branch_counterexample_inverse':['-1/14','1/6'],
 'rational_branch_counterexample_quartic_sides':'75625/3',
 'rational_branch_counterexample_cosine':'-9/5',
 'regular_finite_r_zero_endpoint_formulas_exact':True,
 'scope':'Positive individual maximally entangled C2 generators, including n=0; not their convex hull or a global SDP.',
 'no_numerical_optimization_or_parameter_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
