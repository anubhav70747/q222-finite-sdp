"""Independent full-parent check and rational fixed-gate pencil extraction.

Uses exact SymPy polynomial arithmetic. No SDP solver, sampling, source
envelope sweep, or downloaded matrix data are used.
"""
from pathlib import Path
import json
import sympy as S

HERE=Path(__file__).resolve().parent
a,b,c,e,W,X,Y,Z,g,zeta,rB=S.symbols('a b c e W X Y Z g zeta rB')
d,la,n,m,C,A,D,al,be,tau,p,q=S.symbols('d la n m C A D al be tau p q')
sources=[d,la,n,m,C,A,D,al,be]
source_names=['d','lambda','n','m','C','A','D','alpha','beta']
parameters=sources+[tau,p,q]
B=1-S.Rational(1,2**26); delta=S.Rational(1,2**26)
E=a*W-b*Y; K=a*X-b*Z; J=e*E+c*K; U=e*Z-c*Y
R=b*e*X+a*U; ell=c*Y; V=c*K; Y0=e*Y
El=e*W+c*X; LA=b*X+a*Z
M=b*(c*W-e*X)+a*(c*Y-e*Z)
Da=a*a+b*b; Db=c*c+e*e
Pfull=(la*Da*El**2+n*Da*Db*X**2+m*M**2+C*Db*E**2+A*Da*Db*Y**2
       +(D+d/2)*J**2+(al+d/2)*Da*U**2+(be+d/2)*Db*LA**2-d*Da*Db*Z**2/2)
P=(la*El**2+n*X**2+m*M**2+C*E**2+A*Y**2
   +(D+d/2)*J**2+(al+d/2)*U**2+(be+d/2)*LA**2-d*Z**2/2)

def reduce_circle(expr,kind='alice'):
    expr=S.Poly(S.expand(expr),g,zeta,rB)
    out=0
    for (ig,iz,ir),coef in expr.terms():
        out+=coef*g**(ig%2)*(b*(B-b))**(ig//2)*zeta**(iz%2)*(1-b)**(iz//2)*rB**(ir%2)*(c*c-delta*delta)**(ir//2)
    out=S.Poly(S.expand(out),a,c)
    result=0
    for (ia,ic),coef in out.terms():
        result+=coef*a**(ia%2)*(1-b*b)**(ia//2)*c**(ic%2)*(1-e*e)**(ic//2)
    return S.expand(result)

assert reduce_circle(Pfull-P)==0
assert reduce_circle(LA**2-Z**2-X**2+K**2)==0
Stau=S.Matrix([[C-tau,tau],[tau,D-tau]])
pair=S.Matrix([e*E,J])
retained=la*El**2+be*LA**2+n*X**2+m*M**2+C*(c*E)**2+(pair.T*Stau*pair)[0]+d*J**2/2
rho=tau-d/2
G=d*R**2/2+al*U**2+A*ell**2+d*(b*ell)**2/2+d*(c*X)**2/2+rho*V**2+A*Y0**2+d*b*V*Y0
assert reduce_circle(P-retained-G)==0
assert reduce_circle(b*V*Y0-a*R*ell+U*ell+b*b*ell**2)==0

# Expanded form of (3), using unrestricted independent p,q in place of
# their later nonlinear source-only witness assignments.
Gallocated=(d*R**2/2+p*R*(a*b*ell)+q*(a*b*ell)**2/2
            +al*U**2-p*U*(b*ell)+(d/2-p)*(b*ell)**2+A*ell**2
            +(p-q)*(b*zeta*ell)**2+q*(b*(1-b)*ell)**2/2+d*(c*X)**2/2
            +rho*V**2+A*Y0**2+(d*b-p*b*b)*V*Y0)
assert reduce_circle(Gallocated-G)==0
Hplus=S.Matrix([[rho,(d*B-p*B*B)/2],[(d*B-p*B*B)/2,A]])
H0=S.diag(rho,A)
Qstar=S.Matrix([[2*rho/B**2,d/(2*B)],[d/(2*B),2*A/B**2]])
interpolation=b*b/B**2*Hplus+(B-b)**2/B**2*H0+b*(B-b)*Qstar
assert interpolation.applyfunc(S.expand)==S.Matrix([[rho,(d*b-p*b*b)/2],[(d*b-p*b*b)/2,A]])

HB=S.Matrix([El,LA,X,M,c*E,e*E,J,R,a*b*ell,U,b*ell,ell,
             b*zeta*ell,b*(1-b)*ell,c*X,Y0,b*V,b*Y0,g*V,g*Y0])
Q20=S.zeros(20)
def addblock(Q,indices,mat):
    for i,ii in enumerate(indices):
        for j,jj in enumerate(indices):Q[ii,jj]+=mat[i,j]
def outer(Q,entries,weight):
    for i,x in entries.items():
        for j,y in entries.items():Q[i,j]+=weight*x*y
for index,weight in [(0,la),(1,be),(2,n),(3,m),(4,C),(11,A),(12,p-q),(13,q/2),(14,d/2)]:
    Q20[index,index]+=weight
addblock(Q20,[5,6],Stau)
Q20[6,6]+=d/2
addblock(Q20,[7,8],S.Matrix([[d/2,p/2],[p/2,q/2]]))
addblock(Q20,[9,10],S.Matrix([[al,-p/2],[-p/2,d/2-p]]))
addblock(Q20,[16,17],Hplus/B**2)
outer(Q20,{5:-1,6:1,16:-1/B},rho)
outer(Q20,{15:1,17:-1/B},A)
addblock(Q20,[18,19],Qstar)
Q20=Q20.applyfunc(S.expand)
assert Q20==Q20.T
alice_difference=S.expand((HB.T*Q20*HB)[0]-Pfull)
for parameter in parameters:
    assert reduce_circle(alice_difference.coeff(parameter))==0,str(parameter)
assert S.expand(alice_difference.subs(dict.fromkeys(parameters,0)))==0

tau0=d/(2*delta**2)
Psi=S.Matrix([El,LA,X,M,Y,U,J,c*E,e*E,rB*K])
Q10=S.zeros(10)
for index,weight in [(0,la),(1,be),(2,n+d/2),(3,m),(4,A),(5,al+d/2),(6,d/2),(7,C),(9,tau0)]:
    Q10[index,index]+=weight
addblock(Q10,[8,6],S.Matrix([[C-tau0,tau0],[tau0,D-tau0]]))
Q10=Q10.applyfunc(S.expand)
bob_difference=S.expand((Psi.T*Q10*Psi)[0]-Pfull)
for parameter in sources:
    assert reduce_circle(bob_difference.coeff(parameter))==0,str(parameter)
assert reduce_circle(C*(c*E)**2+(pair.T*S.Matrix([[C-tau0,tau0],[tau0,D-tau0]])*pair)[0]
                     +tau0*(c*c-delta**2)*K**2-C*E**2-D*J**2+d*K**2/2)==0

def export(Q,pars,names,features,filename,domain):
    entries=[]
    for k,parameter in enumerate(pars):
        for i in range(Q.rows):
            for j in range(i,Q.cols):
                coef=S.expand(Q[i,j]).coeff(parameter)
                if coef:
                    assert coef.is_Rational
                    entries.append([k,i,j,str(coef)])
    reconstructed=S.zeros(Q.rows)
    for k,i,j,coef in entries:
        reconstructed[i,j]+=S.Rational(coef)*pars[k]
        if i!=j:reconstructed[j,i]+=S.Rational(coef)*pars[k]
    assert reconstructed==Q
    result={'order':Q.rows,'source_order':source_names,'parameter_order':names,
            'auxiliary_order':names[9:],'source_nonnegativity_external':True,
            'canonical_order_with_nine_source_sign_blocks':Q.rows+9,
            'indexing':'zero-based upper-triangular COO; mirror off-diagonal once',
            'entry_order':['parameter','row','column','rational_coefficient'],
            'entries':entries,'functions':[S.sstr(x) for x in features],
            'B':str(B),'delta':str(delta),'testing_domain':domain,
            'target':'The unmultiplied compact homogeneous parent F_h/4',
            'scope':'An exact certificate cone on this fixed testing chart; global canonical adequacy is not asserted.'}
    (HERE/filename).write_text(json.dumps(result,indent=2)+'\n')
    return len(entries)

na=export(Q20,parameters,source_names+['tau','p','q'],HB,'alice_order20.coo.json',
          ['a^2+b^2=1','c^2+e^2=1','0<=b<=B','g>=0','g^2=b(B-b)','zeta>=0','zeta^2=1-b'])
nb=export(Q10,sources,source_names,Psi,'bob_order10.coo.json',
          ['a^2+b^2=1','c^2+e^2=1','b>=0','c^2>=delta^2','rB>=0','rB^2=c^2-delta^2'])
assert (na,nb)==(36,15)

# Exact source witness factorization behind (3), including rank drops.
h=S.symbols('h')
theta=h/(1+h)
ps=d*theta; qs=d*theta**2
assert S.factor((d*qs-ps**2)/4)==0
second_det=S.factor(al*(d/2-ps)-ps**2/4)
assert S.factor(second_det-d*(2*al*(1-h*h)-d*h*h)/(4*(1+h)**2))==0
result={'full_homogeneous_parent_normalization':True,'full_parent_residual_exact':True,
        'allocation_identity_2':True,'allocation_identity_3_with_unrestricted_p_q':True,
        'interpolation_identity_4':True,'alice_all_12_coefficient_identities':True,
        'alice_upper_triangular_records':na,'bob_identity_5':True,
        'bob_all_9_coefficient_identities':True,'bob_upper_triangular_records':nb,
        'source_witness_rank_one_factorizations_exact':True,
        'testing_endpoints':'Identities hold on the closed algebraic circles; no inversion of a,b,c,e or a state coordinate occurs. Alice includes b=0,B; Bob includes |c|=delta and c=+/-1,e=0.',
        'scope':'The checker verifies the full-parent identities. Analytic source bounds are proved in the PDF. The contact cap and global source coverage remain unresolved.'}
(HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
