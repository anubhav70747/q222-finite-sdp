"""Exact two-parameter C2 hull and dual coefficient identities; no solver."""
import json
from pathlib import Path
import sympy as s

r,x=s.symbols('r x',positive=True)
z=x/r
k=b=s.Integer(1)
t=1-r
h=b-k*t
Ac=1+z*r
g=(1+z)*b-(1-r+k*z*t)
j=1+z*r-t
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+k*z)*b
       +z*z*(1-k)*r*r+2*z*(1+k*z)*r+1+k*z)
N=(1+z)*(b+r)*(b+k*(1-t))-k*(1+k*z)*(1-r)*t
native=[1,0,((b+r)*Ac*(Ac+k*z*t)-r*Delta)/(2*r*Delta),
        (1+k*z)*(1-r)*Ac*(Ac+k*z*t)/(2*g*Delta),
        Ac*N/(2*k*t*Delta),k*Ac*(b+r)*(Ac+k*z*t)/(2*h*Delta),
        (1+z)*(b+r)*(b-k*z*r)*(Ac+k*z*t)/(2*k*j*Delta),
        z*z*(b+r)*N/(2*(1-z*h)*Delta),0]
S=r+x+2*r*x
surface=[1,0,(1-r-x-3*r*x)/(4*S),(1-r)*(1+x)/(4*S),
         r*r*(1+x)/((1-r)*S),(1+r)*(1+x)/(4*S),
         (1+r)*(1-x)/(4*S),x*x*(1+r)/((1-x)*S),0]
assert all(s.cancel(a-b)==0 for a,b in zip(native,surface))
d,la,n,m,C,A,D,al,be=surface
assert s.cancel(d-(4*A-2*m-2*D))==0
assert s.cancel(n-(m+D-A))==0
assert s.cancel(m*C-(A-m)**2)==0
assert s.cancel(D*al-(A-D)**2)==0
assert s.cancel((A-m)/(A+m)-r)==0
assert s.cancel((A-D)/(A+D)-x)==0

AA,mm,DD,CC,aa=s.symbols('A m D C alpha',real=True)
rr=(AA-mm)/(AA+mm)
xx=(AA-DD)/(AA+DD)
assert s.factor(1-rr-xx-3*rr*xx)==4*AA*(mm+DD-AA)/((AA+mm)*(AA+DD))

e=s.symbols('e',positive=True)
Cpath=[s.cancel(s.sympify(v).subs({r:1-e,x:e*e/4})) for v in surface]
Caxis=[s.limit(s.cancel(v/Cpath[4]),e,0,dir='+') for v in Cpath]
apath=[s.cancel(s.sympify(v).subs({x:1-e,r:e*e/4})) for v in surface]
aaxis=[s.limit(s.cancel(v/apath[7]),e,0,dir='+') for v in apath]
assert Caxis==[0,0,0,0,1,0,0,0,0]
assert aaxis==[0,0,0,0,0,0,0,1,0]
assert s.expand((1-r-x-3*r*x).subs({r:1-e,x:e*e/4}))==e-e*e+3*e**3/4

# Full dual polynomial identity, not a fit at sample behaviors.
lA,lm,lD,lC,laa,xi,eta,ss,tt=s.symbols('lA lm lD lC lalpha xi eta s t')
U=s.Matrix([[lA+lm-2*eta-tt,xi],[xi,lC]])
V=s.Matrix([[lA+lD-2*xi-ss,eta],[eta,laa]])
M=s.Matrix([[mm,AA-mm],[AA-mm,CC]])
J=s.Matrix([[DD,AA-DD],[AA-DD,aa]])
q3=2*xi+2*eta+ss+tt-lA
rhs=s.trace(U*M)+s.trace(V*J)+ss*(AA-mm)+tt*(AA-DD)+q3*(mm+DD-AA)
assert s.expand(rhs-(lA*AA+lm*mm+lD*DD+lC*CC+laa*aa))==0

HH,pn,pm,pC,pA,pD,pal=s.symbols('H pn pm pC pA pD palpha')
Bell=(4*AA-2*mm-2*DD)*HH+4*((mm+DD-AA)*pn+mm*pm+CC*pC+AA*pA+DD*pD+aa*pal)
behavior={lA:4*(HH-pn+pA),lm:-2*HH+4*(pn+pm),lD:-2*HH+4*(pn+pD),lC:4*pC,laa:4*pal}
assert s.expand(Bell-(lA*AA+lm*mm+lD*DD+lC*CC+laa*aa).subs(behavior))==0

strict={AA:1,mm:s.Rational(3,4),DD:s.Rational(3,4),CC:1,aa:1}
assert M.subs(strict).det()==s.Rational(11,16)
assert J.subs(strict).det()==s.Rational(11,16)
data={
 'source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
 'native_surface_specialization_exact':True,
 'linear_relations_exact':True,'both_surface_determinants_zero':True,
 'inverse_contact_parameters_exact':True,
 'physicality_domain_identity':'1-r-x-3*r*x = 4*A*(m+D-A)/((A+m)*(A+D))',
 'C_axis_limit':list(map(str,Caxis)),
 'alpha_axis_limit':list(map(str,aaxis)),
 'axis_path_native_domain_gap':'e-e^2+3*e^3/4 > 0 for 0<e<1',
 'strict_source_point':{'d':'1','n':'1/2','m':'3/4','C':'1','A':'1','D':'3/4','alpha':'1','lambda':'0','beta':'0'},
 'strict_matrix_determinants':['11/16','11/16'],
 'strict_scalar_values':['1/4','1/4','1/2'],
 'dual_identity_exact':True,'behavior_pairing_coefficients_exact':True,
 'primal_total_PSD_order':7,'dual_total_PSD_order':7,'dual_auxiliary_count':4,
 'no_solver_or_numerical_sweep':True,
 'scope':'Closed hull of the a=b_source=1, t_source=1-r native C2 surface; not full C2 or full Q222.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8', newline='\n')
print(json.dumps(data,indent=2))
