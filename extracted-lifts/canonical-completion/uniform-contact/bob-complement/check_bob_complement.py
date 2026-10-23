"""Exact Bob-complement identity, rational pencil, and R1 admission.

No optimizer, source sweep, sampling, or numerical matrix test is used.
The supplied published reference is checked before exporting the new
threshold-specific source-affine matrix. Requires Python3 and SymPy.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
a,b,c,e,W,X,Y,Z,rB=s.symbols('a b c e W X Y Z rB')
source=s.symbols('d lambda n m C A D alpha beta')
d,la,n,m,C,A,D,al,be=source
names=['d','lambda','n','m','C','A','D','alpha','beta']
E=a*W-b*Y;K=a*X-b*Z;J=e*E+c*K;U=e*Z-c*Y
El=e*W+c*X;LA=b*X+a*Z
M=b*(c*W-e*X)+a*(c*Y-e*Z)
features=s.Matrix([El,LA,X,M,Y,U,J,c*E,e*E,rB*K])
Da=a*a+b*b;Db=c*c+e*e
parent=(la*Da*El**2+n*Da*Db*X**2+m*M**2+C*Db*E**2+A*Da*Db*Y**2
       +(D+d/2)*J**2+(al+d/2)*Da*U**2+(be+d/2)*Db*LA**2-d*Da*Db*Z**2/2)

def matrix(q):
    tau=d/(2*q)
    Q=s.diag(la,be,n+d/2,m,A,al+d/2,D+d/2-tau,C,C-tau,tau)
    Q[6,8]=Q[8,6]=tau
    return Q

def reduce_circle(expr,q):
    out=0
    for (power,),coef in s.Poly(s.expand(expr),rB).terms():
        out+=coef*rB**(power%2)*(c*c-q)**(power//2)
    ans=0
    for (ia,ic),coef in s.Poly(s.expand(out),a,c).terms():
        ans+=coef*a**(ia%2)*(1-b*b)**(ia//2)*c**(ic%2)*(1-e*e)**(ic//2)
    return s.expand(ans)

refpath=HERE/'published_bob_order10_reference.coo.json'
ref=json.loads(refpath.read_text())
old=s.zeros(10)
for k,i,j,coef in ref['entries']:
    old[i,j]+=s.Rational(coef)*source[k]
    if i!=j:old[j,i]+=s.Rational(coef)*source[k]
qold=s.Rational(ref['delta'])**2
assert old==matrix(qold)
assert len(ref['entries'])==15
local={str(x):x for x in [a,b,c,e,W,X,Y,Z,rB]}
assert all(s.expand(s.sympify(x,locals=local)-f)==0 for x,f in zip(ref['functions'],features))

q=s.Rational(1,257);Q=matrix(q)
difference=s.expand((features.T*Q*features)[0]-parent)
assert all(reduce_circle(difference.coeff(par),q)==0 for par in source)
assert s.expand(difference.subs(dict.fromkeys(source,0)))==0
assert reduce_circle(LA**2-Z**2-X**2+K**2,q)==0
tau=d/(2*q)
S=s.Matrix([[C-tau,tau],[tau,D-tau]])
pair=s.Matrix([e*E,J])
allocation=C*(c*E)**2+(pair.T*S*pair)[0]+tau*(c*c-q)*K*K-C*E*E-D*J*J+d*K*K/2
assert reduce_circle(allocation,q)==0

entries=[]
for k,par in enumerate(source):
    for i in range(10):
        for j in range(i,10):
            coef=s.expand(Q[i,j]).coeff(par)
            if coef:entries.append([k,i,j,str(coef)])
assert len(entries)==15
rebuilt=s.zeros(10)
for k,i,j,coef in entries:
    rebuilt[i,j]+=s.Rational(coef)*source[k]
    if i!=j:rebuilt[j,i]+=s.Rational(coef)*source[k]
assert rebuilt==Q
for par in source[1:]:
    event=Q.diff(par)
    assert all(event[i,j]==0 for i in range(10) for j in range(10) if i!=j)
    assert all(event[i,i]>=0 for i in range(10))
out={'order':10,'source_order':names,'parameter_order':names,'auxiliary_order':[],
 'source_nonnegativity_external':True,'indexing':'Zero-based upper-triangular COO; mirror each off-diagonal once.',
 'entry_order':['parameter','row','column','rational_coefficient'],'entries':entries,
 'functions':list(map(str,features)),'radical_squared':'c^2-1/257',
 'testing_domain':['a^2+b^2=1','c^2+e^2=1','b>=0','c^2>=1/257','rB>=0','rB^2=c^2-1/257'],
 'equivalent_finite_Bob_domain':'|v|>=1/16','target':'Unmultiplied compact canonical parent F_h/4',
 'PSD_admission':'c_i>=0 and [[C-257*d/2,257*d/2],[257*d/2,D-257*d/2]]>=0',
 'scope':'All Alice angles and states on this Bob region, including projective endpoints. Not global canonical adequacy.',
 'joined_order_with_cap':'274 cap +10 Bob +9 shared source signs =293'}
(HERE/'bob_complement_order10.coo.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')

# Native source admission. The source appendix proves C>=k/H and
# D>=(1-delta^2)/(4R0*T*(1-H)); f0<1-delta^2 makes D>=k/(1-H).
delta=F(1,1024);f0=1-delta-delta**3
R0=1+2*delta+2*delta**2+4*delta**3+delta**4*(1+2*delta)
assert f0<1-delta**2
assert f0/(4*R0)>F(1,5)
assert F(1024,5)>F(257,2)
hh,kk,cc,dd,tt=s.symbols('H k C D tau',positive=True)
assert s.simplify(hh/kk+(1-hh)/kk-1/kk)==0
assert s.expand(s.Matrix([[cc-tt,tt],[tt,dd-tt]]).det()-cc*dd+tt*(cc+dd))==0
receipt={'published_reference_SHA256':hashlib.sha256(refpath.read_bytes()).hexdigest(),
 'published_all15_COO_records_match_independent_formula':True,
 'published_all10_functions_match':True,'new_all9_parent_coefficient_identities':True,
 'Bob_allocation_identity':True,'new_all15_rational_COO_records_reconstruct_exactly':True,
 'f0':str(f0),'R0':str(R0),'source_bound_f0_over_4R0':str(f0/(4*R0)),
 'f0_lt_1_minus_delta_squared':True,'f0_over_4R0_gt_1_5':True,
 'eta_strict_lower':'1024/5','tau_native':'257/2','eta_margin_lower':'763/10',
 'literal_zero_faces':'n=0,lambda=0,beta=0 allowed; all eight d=0 event rays PSD; no division by an event coefficient.',
 'projective_endpoints':'b=0,a=+/-1 and e=0,c=+/-1 included algebraically.',
 'remaining_tests':'|u|>9/256 and |v|<1/16, including Alice projective endpoints in that strip.',
 'no_solver_or_numerical_sweep':True}
(HERE/'bob_complement_checks.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(receipt,indent=2))
