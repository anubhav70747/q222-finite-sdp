"""Extract the exact twelve-root Gram from Uniform_C1_Certificate.tex.

No numerical search or guessed Gram frame is used. All source symbols are
defined in source_definitions.json; s=sqrt(u) and w are testing variables.
"""
import json
from pathlib import Path
import sympy as S

HERE=Path(__file__).resolve().parent
s,w=S.symbols('s w',nonnegative=True)
names='b k q sigma p m D beta C a ell h rho dG S0 nu R0 Rx Bx Ay chi lam q0 L0 KA c a0 kappa r'
v=dict(zip(names.split(),S.symbols(names,real=True)))
globals().update(v)
u=s*s; z=q*q; Qb=1-b
n=12
eye=S.eye(n)
def e(i): return eye[:,i]
Nx,Ny,R,GS,At,Tt,d1,d2,d3,d4,Us,Ks=[e(i) for i in range(n)]
def prod(x,y): return (x*y.T+y*x.T)/2
def square(x): return x*x.T

# The source root matrix psi=A phi.  Entries are exact algebraic functions.
phi_names=['X','Y','Z','W','u*Z','w*Z','u*W','u*w*W','u^2*W',
           'u^(3/2)*W','u^(1/2)*W','sqrt(u*w*W*Z)']
A=S.zeros(12)
A[0,0]=1; A[0,4]=h
A[1,1]=1/q; A[1,5]=-b/q
A[2,5]=1/q; A[2,6]=-q
A[3,2]=S.sqrt(dG)*sigma; A[3,7]=-S.sqrt(dG)
A[4,3]=-q*k*k; A[4,8]=q
A[5,3]=-k*k*q; A[5,6]=k*q
A[6,4]=S.sqrt(L0)
A[7,7]=S.sqrt(rho)
A[8,3]=S.sqrt(lam)*k*k*q
A[9,5]=S.sqrt(lam)*k/q
for i in range(6,10): A[i,11]=-S.sqrt(lam)*k
A[10,9]=S.sqrt(q0)*q; A[10,10]=-S.sqrt(q0)*q*k
A[11,6]=S.sqrt(q0)*S.sqrt(k)*q
A[11,3]=-S.sqrt(q0)*S.sqrt(k)*q*k

# Two independent, complete state-normal conductors.
UG=(S.sqrt(dG)*u*GS+S.sqrt(rho)*(k*(d1-d2)-(u-k)*(d2-d3))
    -rho*S.sqrt(sigma)*Tt)/S0
R0wG=(R0*S.sqrt(dG)*q*w*GS+R0*rho*sigma*q*q*R
       -R0*S.sqrt(rho)*(q*w+q*q*S.sqrt(sigma))*(d2-d3)
       +(R0*q0*q*q+R0*S.sqrt(lam)*S.sqrt(rho)*q*w)*Tt)/S0
full_second=R0wG-8*p*q*w*u*UG
JC=Ny+b*R; JL=Ny-Qb*R
base=2*square(Nx)+square(JC)+square(JL)+square(GS)
base+=(KA+c)/2*(square(At)+square(Tt))
base+=sum((square(e(i))/2 for i in range(6,10)),S.zeros(n))
base+=square(Us)+square(Ks)/2

ERR0=z*(8*b*D+8*b*p*sigma+4*b*sigma*(2*C+2*D+1)*(r-1))
EyR=4*z*(u*u*(2*D*h+2*p*sigma+(1+2*m)*(h-Qb))
        +b*(2*D+1)*w*w+sigma*(2*C+1)*(r-1)
        +2*D*(1+sigma*(r-1))+2*p*sigma)
epsC=2/lam+u+k
epsG=2*a0*u+(u+k)/L0
Q=base+(ERR0-lam*k*k/2-c*dG*u*u/S0)*square(R)
Q+=EyR*prod(Ny,R)-4*nu*prod(Nx,UG)+prod(JC,full_second)
Q+=-4*u*q*w*(2*m+2*D+1)*prod(Nx,Ny)
Q+=u*q*w*Rx*prod(Nx,R)+q*w*Bx*prod(Nx,At)
Q+=u*Ay*prod(Ny,At)+chi*u*prod(At,R)
Q+=S.sqrt(q0)*s*prod(Us,R)
Q+=-c/S0*square(S.sqrt(a0)*s*Us+S.sqrt(a0*k)*Ks)
Q+=-c*epsG/S0*square(GS)
Q+=-c*dG*epsC/S0*sum((square(e(i)) for i in range(6,10)),S.zeros(n))
Q=Q.applyfunc(S.expand)

coeffs={}
for i in range(n):
 for j in range(i,n):
  for powers,coef in S.Poly(Q[i,j],s,w).terms():
   coeffs.setdefault(powers,[]).append([i,j,S.sstr(S.factor(coef))])
data={'schema':'Q222.C1.twelve_root_gram.v1','root_order':['Nx','Ny','R','Gamma_star','A_theta','T_theta','d1','d2','d3','d4','U_star','K_star'],
      'testing_parameters':['s=sqrt(u)','w'],'matrix_order':12,
      'convention':'Only i<=j stored; off-diagonal coefficients are matrix entries, not doubled polynomial coefficients.',
      'coefficient_matrices':[{'monomial_exponents':list(k),'entries':coeffs[k]} for k in sorted(coeffs)],
      'source_symbols':list(v)}
(HERE/'C1_PRINCIPAL_GRAM.json').write_text(json.dumps(data,indent=2)+'\n')
Aminus=A.copy()
for i in range(6,10): Aminus[i,11]=0
(HERE/'C1_ROOT_MATRIX.json').write_text(json.dumps({'root_order':data['root_order'],'testing_span':phi_names,
    'source_matrix':[[S.sstr(x) for x in A.row(i)] for i in range(12)],
    'source_matrix_negative_state_sign':[[S.sstr(x) for x in Aminus.row(i)] for i in range(12)],
    'scope':'lambda>0, u,w>=0. Positive state chart uses W,Z>=0; negative state chart uses actual Z<=0 and W>=0, with d_i=B_i and the twelfth column identically zero. All testing axes are included by polynomial specialization.'},indent=2)+'\n')
defs={'x':'k^2','z':'q^2','Q':'1-b','theta':'1-h','H':'h/sigma','e':'k^2*h/sigma',
      'L0':'lambda*sigma','rho':'lambda*k^2/sigma','q0':'lambda*k','dG':'4*D/sigma^2',
      'S0':'dG+rho','a0':'k/sigma','c':'max(0,1-KA)',
      'KA':'4*kappa*H-lambda','a':'p/b','ell':'l/(1-b)',
      'r':'b*z','R0':'-4*(2*(C+D)*(r-1)+2*p+r-1)',
      'nu':'(2*D-2*beta*(1-h)+h)/sigma',
      'lambda':'(8*D*h+4*beta*(1-h)^2+2*h^2)/sigma-4*b*p*sigma*z',
      'Rx':'-4*(2*D*b+(1+2*m)*(h-(1-b))-2*beta*(1-h))',
      'Bx':'8*kappa*e','Ay':'8*(p-l)+4*z*((1+2*m)*(h-(1-b))+2*D*h+2*p*sigma)',
      'chi':'8*(p*b+l*(1-b))-4*h+4*e*h*(2*D+1)-8*D*h+4*z*(4*D*b*h+4*p*b*sigma+(1+2*m)*(h-(1-b))^2+2*beta*(1-h)^2)'}
(HERE/'source_definitions.json').write_text(json.dumps(defs,indent=2)+'\n')
stats={'matrix_order':12,'testing_monomials':[list(k) for k in sorted(coeffs)],
       'stored_upper_entries':sum(len(x) for x in coeffs.values()),
       'source_root_matrix_nonzero_entries':sum(x!=0 for x in A),'base_minimum_eigenvalue_bound':'1/3',
       'base_assumptions':['0<=b<=1','KA+c>=1']}
(HERE/'extraction_summary.json').write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
