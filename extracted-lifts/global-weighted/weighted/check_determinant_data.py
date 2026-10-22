#!/usr/bin/env python3
"""Independent exact verification of all determinant coefficients SG5--SG6.

The determinant is computed by Cauchy--Binet directly from the eight original
affine event/core rows.  No sampling, guessed Gram frame or old pencil is used.
"""
from itertools import combinations, permutations
from pathlib import Path
import json
import sympy as s

u,v=s.symbols('u v')
d,lam,n,m,C,A,D,alpha,beta=s.symbols('d lambda n m C A D alpha beta')
source=(d,lam,n,m,C,A,D,alpha,beta)
variables=(u,v)+source
Du,Dv=1+u*u,1+v*v
R=s.Matrix([[1,v,0,0],[0,1,0,0],[v,-1,u*v,-u],[u,0,-1,0],
            [0,0,1,0],[u,u*v,-1,-v],[0,0,v,-1],[0,1,0,u]])
weights=[Du,Du*Dv,1,Dv,Du*Dv,1,Du,Dv]
kappa=[4*lam,4*n,4*m,4*C,4*A,4*D+2*d,4*alpha+2*d,4*beta+2*d]

def poly(expr):
    return s.Poly(expr,*variables,domain=s.QQ)

detP=poly(0)
for inds in combinations(range(8),4):
    minor=s.expand(R[list(inds),:].det())
    if minor:
        detP += poly(s.prod(kappa[i]*weights[i] for i in inds)*minor**2)
for inds in combinations(range(8),3):
    minor=s.expand(R[list(inds),:3].det())
    if minor:
        detP -= poly(2*d*Du*Dv*s.prod(kappa[i]*weights[i] for i in inds)*minor**2)

t1=d+2*(lam+n+D); t2=d+2*(alpha+m+A)
t3=d+2*(lam+C+D); t4=d+2*(alpha+beta+m)
t5=d+2*(lam+beta+n); t6=d+2*(alpha+C+A)
t7=d+2*(beta+n+m); t8=d+2*(C+A+D)
q={(0,0):4*lam*alpha*t7*t8,
   (0,1):4*(lam*D*t7*t6+alpha*m*t5*t8),
   (0,2):4*m*D*t5*t6,
   (1,0):4*(alpha*A*t7*t3+lam*n*t8*t4),
   (1,2):4*(beta*m*t1*t6+C*D*t5*t2),
   (2,0):4*n*A*t3*t4,
   (2,1):4*(beta*n*t3*t2+C*A*t1*t4),
   (2,2):4*beta*C*t1*t2}

P11=sum((kappa[i]*weights[i].subs({u:1,v:1})*R[i,:].T*R[i,:]
         if hasattr(weights[i],'subs') else kappa[i]*R[i,:].T*R[i,:]
         for i in range(8)),s.zeros(4)).subs({u:1,v:1})
P11[3,3]-=8*d
# A direct 24-term determinant, independent of the exterior-power expansion.
def parity(p):
    return (-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
det11=s.expand(sum(parity(p)*s.prod(P11[i,p[i]] for i in range(4))
                  for p in permutations(range(4))))
h=s.expand(det11/256-sum(q.values()))
q[(1,1)]=h
target=poly(16*Du**2*Dv**2*sum(c*u**(2*i)*v**(2*j) for (i,j),c in q.items()))
assert detP==target

receipt={"all_eight_outer_coefficients_verified":True,
         "central_coefficient_verified":True,
         "determinant_identity_verified_as_polynomial":True,
         "source_variables":[str(x) for x in source],
         "determinant_polynomial_terms":len(detP.terms()),
         "central_coefficient_terms":len(s.Poly(h,*source).terms()),
         "method":"Exact Cauchy--Binet from the eight affine rows; central coefficient independently obtained from the 24-term determinant at u=v=1.",
         "old_order669_recomputed":False}
data={"variables":[str(x) for x in source],
      "coefficients":[{"ij":[i,j],"terms":[[list(e),int(c.p),int(c.q)]
          for e,c in s.Poly(q[i,j],*source).terms()]}
          for i,j in sorted(q)],
      "rule":"det P_aff = 16(1+u^2)^2(1+v^2)^2 sum q_ij u^(2i)v^(2j); h=q_11."}
HERE=Path(__file__).resolve().parent
(HERE/'DETERMINANT_SOURCE_COEFFICIENTS.json').write_text(json.dumps(data,indent=2)+'\n')
(HERE/'determinant_check_results.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
