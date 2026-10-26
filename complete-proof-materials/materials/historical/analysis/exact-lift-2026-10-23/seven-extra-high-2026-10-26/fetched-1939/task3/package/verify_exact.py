#!/usr/bin/env python3
"""Self-contained exact proof replay, corrected C2 seed.

NO SDP solver, numerical search, samples, external files or saved intermediates.
It proves the full polynomial syzygy contact-jet capacity, a source-only
rational 32-packet two-budget exclusion when the old tau=2 top degree packet
is frozen, and shows that releasing those rows removes that exclusion.
This is NOT a counterexample to K364, nor a universal K364 certificate.
"""
import json
from pathlib import Path
import sympy as sp
Q=sp.Rational
u,v,W,X,Y,Z=sp.symbols('u v W X Y Z')
ct={u:1,v:1}
n,m,C,A,D,al=Q(1,8),Q(3,8),Q(1,6),Q(5,8),Q(3,8),Q(1,6)
row=[]; k2=[]
for f in (1,u):
    for g in (1,v): row.append((0,f*g,0,0));k2.append(4*n)
row.append((v,-1,u*v,-u));k2.append(4*m)
for g in (1,v): row.append((u*g,0,-g,0));k2.append(4*C)
for f in (1,u):
    for g in (1,v): row.append((0,0,f*g,0));k2.append(4*A)
row.append((u,u*v,-1,-v));k2.append(4*D+2)
for f in (1,u): row.append((0,0,f*v,-f));k2.append(4*al+2)
for g in (1,v): row.append((0,g,0,u*g));k2.append(Q(2))
assert len(row)==16
R=sp.Matrix(row);R=R.col_join(v*R);k2=k2+k2
psi=sp.Matrix([1,-Q(1,4),Q(1,4),1]); KK=sp.diag(*k2);Wgt=sp.diag(*[1/k for k in k2])
P=sum(k2[i]*(R[i,:]*sp.Matrix([W,X,Y,Z]))[0]**2 for i in range(32))-2*(1+u*u)*(1+v*v)**2*Z**2
assert all(sp.simplify(expr.subs({**ct,W:1,X:-Q(1,4),Y:Q(1,4),Z:1}))==0 for expr in [P]+[sp.diff(P,z) for z in (W,X,Y,Z,u,v)])
Ppsi=sp.expand(P.subs({W:1,X:-Q(1,4),Y:Q(1,4),Z:1}))
assert sp.factor(Ppsi.subs(v,1))==Q(289,24)*(u-1)**2
assert sp.factor(Ppsi.subs(u,0).subs(v,0))>=0
print('PASS: exact complete four-state parent, all six contact jets')

# Normal syzygies in the user's original positive-packet basis (32 rows).
S=sp.zeros(32,29); j=0
for i in (1,2,3,14,15):S[i,j]=1;S[0,j]=-R[i,1];j+=1
for i in (8,9,10,12,13):S[i,j]=1;S[7,j]=-R[i,2];j+=1
for entry in ({6:1,5:-v},{11:1,5:-1,0:-u*v},{4:u,5:-v,0:u,7:-v,9:-u*v}):
    for i,z in entry.items():S[i,j]=z
    j+=1
for i in range(16):S[16+i,j]=1;S[i,j]=-v;j+=1
assert j==29 and R[:,:3].T*S==sp.zeros(3,29)
high=sp.zeros(32,1)
for i,z in {27:v*v,20:-u*v,31:-u*v*v,16:-u*v,24:v,26:u*v}.items():high[i]=z
assert R[:,:3].T*high==sp.zeros(3,1)

def interp(b):
    eta=[b[i]+b[16+i] for i in range(16)]
    inds=(1,2,3,14,15,8,9,10,12,13,6,11,4)
    return sp.Matrix([eta[i] if i!=9 else eta[9]+eta[4] for i in inds]+list(b[16:]))

b=-sp.Matrix([k2[i]*(R[i,:]*psi)[0].subs(ct) for i in range(32)])
t0=S*interp(b)+2*(high-S*interp(high.subs(ct)))
assert t0.subs(ct)==b
# Check the *complete four-state* polynomial identity for the original
# contact-interpolated source-only syzygy, not only its state restriction.
normal_roots=[(k2[i]*(R[i,:]*sp.Matrix([W,X,Y,0]))[0]
               +(k2[i]*R[i,3]+t0[i])*Z) for i in range(32)]
scalar0=-2*(1+u*u)*(1+v*v)**2-2*sum(R[i,3]*t0[i] for i in range(32))-sum(t0[i]**2/k2[i] for i in range(32))
assert sp.expand(P-sum(normal_roots[i]**2/k2[i] for i in range(32))-scalar0*Z**2)==0
G0=sp.Matrix([sp.expand(k2[i]*(R[i,:]*psi)[0]+t0[i]) for i in range(32)])
gu=G0.diff(u).subs(ct)
assert sp.factor((gu.T*Wgt*gu)[0])==Q(1523,40)

# Entire polynomial syzygy module: t_i=sum_{p=0}^1 sum_{q=0}^2 z[i,p,q]u^p v^q.
mon=[(i,p,q) for i in range(32) for p in range(2) for q in range(3)]
Aeq=[]
for state in range(3):
    for pu in range(3):
        for pv in range(5):
            Aeq.append([sp.Poly(R[i,state],u,v).coeff_monomial(u**(pu-p)*v**(pv-q)) if pu>=p and pv>=q else 0 for i,p,q in mon])
Aeq=sp.Matrix(Aeq)
E=sp.Matrix([[int(i==ii) for ii,p,q in mon] for i in range(32)])
U=sp.Matrix([[p*int(i==ii) for ii,p,q in mon] for i in range(32)])
Hu=sp.Matrix([[int(ii==i and p==1 and q==2) for ii,p,q in mon] for i in range(32)])
Hv=sp.Matrix([[int(ii==i and p==0 and q==2) for ii,p,q in mon] for i in range(32)])
Ou=sp.Matrix([[int(ii==i and p==1 and q==0) for ii,p,q in mon] for i in range(32)])
Oo=sp.Matrix([[int(ii==i and p==0 and q==0) for ii,p,q in mon] for i in range(32)])
assert Aeq.shape==(45,192)
# Necessary and sufficient constraints on the u-jet of contact-zero syzygies.
null=(Aeq.col_join(E).col_join(U)).T.nullspace()
Lbig=sp.Matrix([[col[77+i] for i in range(32)] for col in null])
L=sp.Matrix.hstack(*Lbig.T.columnspace()).T
assert L.rank()==6 and L.shape==(6,32)
Gu=KK*L.T*(L*KK*L.T).inv()*(L*gu)
assert sp.factor((Gu.T*Wgt*Gu)[0])==Q(289,24)
print('PASS: complete-degree u-jet image has six annihilators; exact optimal curvature 289/24')

# Freeze the old tau=2 high Bob-degree allocation. This is an ansatz-only
# restriction, NOT a restriction of general K364 Grams.
Csys=Aeq.col_join(E).col_join(U).col_join(Hu).col_join(Hv)
rhs=sp.zeros(45,1).col_join(sp.zeros(32,1)).col_join(Gu-gu).col_join(sp.zeros(64,1))
sol,par=Csys.gauss_jordan_solve(rhs)
base=sol.subs({p:0 for p in par})
N=sp.Matrix([[sp.diff(sol[i],p) for p in par] for i in range(sol.rows)])
assert len(par)==38
T=Ou.col_join(Oo)
BB=T*N
# Along v=0 and psi=psi_c, the original full parent has leading coefficient
# 409/48 for u^2; at the origin with psi=(0,0,0,1), it equals 2/3.
f_u=sp.Matrix([sp.Poly(sp.expand(k2[i]*(R[i,:]*psi)[0]+t0[i]),u,v).coeff_monomial(u) for i in range(32)])
f_o=sp.Matrix([k2[i]*R[i,3].subs({u:0,v:0})+t0[i].subs({u:0,v:0}) for i in range(32)])
d=f_u.col_join(f_o)+T*base
WW=sp.diag(*([1/k for k in k2]*2))
Gram=BB.T*WW*BB
assert Gram.det()!=0
minnorm=sp.factor((d.T*WW*d)[0]-(d.T*WW*BB*Gram.inv()*BB.T*WW*d)[0])
assert minnorm==Q(57818093029,6288328800)
cap=Q(409,48)+Q(2,3)
margin=sp.factor(cap-minnorm)
assert margin==-Q(44072179,6288328800)
print('PASS: FROZEN-top packet has negative exact combined budget',margin)

# The same necessary budget is NOT obstructed when the degree-v^2 rows are
# allowed to move: this rules out treating the above as a K364 separator.
relaxed={}
for label,extra,need_min in [('all',[],Q(1019,144)),('only_uv2',[Hu],Q(233201,29952))]:
    mat=Aeq.col_join(E).col_join(U)
    y=sp.zeros(45,1).col_join(sp.zeros(32,1)).col_join(Gu-gu)
    for extra_mat in extra:mat=mat.col_join(extra_mat);y=y.col_join(sp.zeros(32,1))
    sx,fp=mat.gauss_jordan_solve(y)
    xbase=sx.subs({p:0 for p in fp})
    NN=sp.Matrix([[sp.diff(sx[i],p) for p in fp] for i in range(sx.rows)])
    Bbig=T*NN
    drel=f_u.col_join(f_o)+T*xbase
    # Dual projection uses only 9 or 15 exact rational pivots.
    Uperp=sp.Matrix.hstack(*Bbig.T.nullspace())
    invWW=sp.diag(*([z for z in k2]*2))
    Gperp=Uperp.T*invWW*Uperp
    norm=sp.factor((drel.T*Uperp*Gperp.inv()*Uperp.T*drel)[0])
    assert norm==need_min
    relaxed[label]=str(cap-norm)
assert Q(relaxed['all'])>0 and Q(relaxed['only_uv2'])>0
print('PASS: releasing the old high-degree allocation removes this two-budget obstruction')

record=dict(source_native=['1','1','1','3/4','1/4'],weights=['1/8','3/8','1/6','5/8','3/8','1/6'],
     initial_u_jet_energy='1523/40',sharp_u_jet_energy='289/24',
     frozen_high_degree_two_budget_max=str(margin),
     two_budget_unfrozen_margins=relaxed,
     scope='Ansatz-only frozen high-Bob packet no-go. No obstruction to K364; H2 open.')
Path(__file__).with_name('EXACT_RESULTS.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS: ALL EXACT CLAIMS; universal K364 admission NOT established')
