#!/usr/bin/env python3
"""Generate and verify ONE exact rational finite-c native C2 Gram, using proved symbolic recipe.

No numerical optimization and no sweeps. A source-only dyadic gauge terminates
because of the strictly positive analytical Schur limit. An input epsilon must
be small enough for the returned Gram to be positive definite on 12+1 roots.
Examples of supported usage:
python construct_exact_witness.py --a 1 --t 3/2 --eta 1/4 --rho 1 --epsilon 1/1000 --out gram.json

The example above is only a smoke test: the theorem is proved symbolically
in FINITE_C_THEOREM_AND_PROOF.md and verify_finite_c.py for every strict face.
"""
import argparse,json,sys
import sympy as s
import verify_finite_c as v

def rational(x):return s.Rational(x)
def positive_definite(M):
 """Exact rational Sylvester/Schur pivots, never floating-point eigenvalues."""
 A=s.Matrix(M);n=A.rows;piv=[]
 for k in range(n):
  z=s.cancel(A[k,k]);piv.append(z)
  if z<=0:return False,piv
  for i in range(k+1,n):
   for j in range(i,n):
    A[i,j]=s.cancel(A[i,j]-A[i,k]*A[k,j]/z)
  for i in range(k+1,n):
   for j in range(k+1,i):A[i,j]=A[j,i]
 return True,piv

def execute(args):
 sub={v.a:rational(args.a),v.t:rational(args.t),v.eta:rational(args.eta),v.rho:rational(args.rho)}
 aa,tt,ee,rr=[sub[i] for i in (v.a,v.t,v.eta,v.rho)]
 assert aa>0 and tt>0 and 0<rr<tt and tt-rr<1 and 0<ee<1 and aa**2*tt+ee>1
 epv=rational(args.epsilon);assert epv>0
 # Select the first permitted dyadic pure-source gauge with positive two-by-two Schur.
 maxkap=s.cancel((4*v.a*v.p2*v.q2/(v.q2+v.a**2*v.p2)).subs(sub))
 Be=v.Btail[:5,:].subs(sub)
 T=s.diag(v.TX.subs(sub),v.TY.subs(sub))
 for gauge_index in range(2,args.max_gauge_index):
  kx=maxkap/s.Integer(2)**gauge_index
  Ac=v.Gk[:5,:5].subs(sub).subs(v.kap,kx)
  # A_c PD is guaranteed from 0<kx<kmax by the source-wide endpoint LDL.
  S2=(T-Be.T*Ac.inv()*Be).applyfunc(s.cancel)
  if S2[0,0]>0 and s.factor(S2.det())>0:break
 else:raise RuntimeError('Gauge cap exhausted (increase --max-gauge-index); theorem guarantees eventual termination.')
 # A fixed 66-pivot rational right inverse; no original binary intermediates.
 nr={key:s.cancel(val.subs(sub).subs(v.eps,epv)) for key,val in v.Ws.items()}
 fraw=v.fv
 fexact=s.Matrix([s.cancel(expr.subs(sub).subs(v.eps,epv).subs(nr)/4) for expr in fraw])
 f0=v.f0.subs(sub);f1=v.f1.subs(sub)
 q=s.Matrix([s.cancel((fexact[i]-f0[i]-epv*f1[i])/epv**2) for i in range(78)])
 hz=v.Rinv*q[v.pr2,:]
 G2=s.zeros(79,1)
 for j,col in enumerate(v.pc2):G2[col]=hz[j]
 G=(v.H0.subs(sub).subs(v.kap,kx)+epv*v.H1.subs(sub)+epv**2*G2).applyfunc(s.cancel)
 assert all(s.cancel(x)==0 for x in v.M*G-fexact), 'Exact full-parent coefficient identity failed'
 even=s.zeros(12)
 for z,(i,j) in enumerate(v.DIFFCOLS):
  even[i,j]=G[z];even[j,i]=G[z]
 pd,piv=positive_definite(even)
 odd=G[78]
 status='PSD' if pd and odd>0 else 'EXACT_IDENTITY_BUT_EPSILON_NOT_SMALL_ENOUGH_FOR_SELECTED_GAUGE'
 if args.out:
  o={'source':{str(k):str(val) for k,val in sub.items()},'epsilon':str(epv),
     'gamma':str(kx),'gamma_index':gauge_index,'status':status,
     'odd_Gram_scalar':str(odd),'exact_12x12_Gram_upper':[
       [i,j,str(even[i,j])] for i in range(12) for j in range(i,12) if even[i,j]],
     'LDL_pivots_until_first_failure':[str(q) for q in piv],
     'source_weights':{str(k):str(z) for k,z in nr.items()},
     'coefficient_identity':'EXACTLY_VERIFIED_IN_SYMPY'}
  with open(args.out,'w') as f:json.dump(o,f,indent=2)
 print('Exact polynomial identity: VERIFIED; source-only gamma index:',gauge_index, 'status:',status)
 if status!='PSD':return 2
 return 0
if __name__=='__main__':
 p=argparse.ArgumentParser()
 for k in ['a','t','eta','rho','epsilon']:p.add_argument('--'+k,required=True)
 p.add_argument('--max-gauge-index',type=int,default=100)
 p.add_argument('--out',default='')
 sys.exit(execute(p.parse_args()))
