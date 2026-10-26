#!/usr/bin/env python3
"""Exact, self-contained C2 K364 polynomial normal-transfer audit.

Checks the 16-row original positive packet; the (1+v^2) two-copy
32-row parent; 29 explicit polynomial syzygies of its full W/X/Y
state matrix with (u,v) degrees at most (1,2); and contact-kernel
interpolation by pure-source scalings. No solver, samples, network,
or absent runtime files. Does NOT prove the remaining scalar SOS.
"""
import sympy as s
import json
from pathlib import Path
u,v,W,X,Y,Z=s.symbols('u v W X Y Z', real=True)
n,m,C,A,D,al=s.symbols('n m C A D alpha', positive=True)
Du,Dv=1+u*u,1+v*v
rows=[]
def put(name,row,k): rows.append((name,tuple(map(s.expand,row)),k))
for f,fs in [(1,'1'),(u,'u')]:
    for g,gs in [(1,'1'),(v,'v')]:
        put('n_'+fs+gs,(0,f*g,0,0),2*s.sqrt(n))
put('m',(v,-1,u*v,-u),2*s.sqrt(m))
for g,gs in [(1,'1'),(v,'v')]:
    put('C_'+gs,(u*g,0,-g,0),2*s.sqrt(C))
for f,fs in [(1,'1'),(u,'u')]:
    for g,gs in [(1,'1'),(v,'v')]:
        put('A_'+fs+gs,(0,0,f*g,0),2*s.sqrt(A))
put('D',(u,u*v,-1,-v),s.sqrt(4*D+2))
for f,fs in [(1,'1'),(u,'u')]:
    put('alpha_'+fs,(0,0,f*v,-f),s.sqrt(4*al+2))
for g,gs in [(1,'1'),(v,'v')]:
    put('core_'+gs,(0,g,0,u*g),s.sqrt(2))
assert len(rows)==16
kappa=[item[2] for item in rows]
base=s.Matrix([list(z[:3]) for _,z,_ in rows])
fullrows=s.Matrix([list(z) for _,z,_ in rows])
F_pos=sum(k*k*sum(zz*xx for zz,xx in zip(z,(W,X,Y,Z)))**2 for _,z,k in rows)
F_expected=(4*n*Du*Dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2
      +4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2
      +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2
      +(4*al+2)*Du*(v*Y-Z)**2+2*Dv*(X+u*Z)**2)
assert s.expand(F_pos-F_expected)==0
assert s.expand(F_pos-2*Du*Dv*Z**2-
       (F_expected-2*Du*Dv*Z**2))==0
print('PASS: full four-state positive-packet identity, with original negative Z square')

# Baseline simplified row system indexed by X0=0, m=4, C0=5, Y0=7, Yu=9.
X0, m_idx,C0,Y0,Yu=0,4,5,7,9
cols=[]
def syz(name,co):
    vct=s.zeros(32,1)
    for i,ex in co.items():vct[i]=s.expand(ex)
    cols.append((name,vct))
for i in [1,2,3,14,15]:
    factor=base[i,1]
    syz('X_'+str(i),{i:1,X0:-factor})
for i in [8,9,10,12,13]:
    factor=base[i,2]
    syz('Y_'+str(i),{i:1,Y0:-factor})
syz('Cv',{6:1,C0:-v})
syz('D',{11:1,C0:-1,X0:-u*v})
syz('m',{m_idx:u,C0:-v,X0:u,Y0:-v,Yu:-u*v})
for i in range(16):
    syz('shift_'+str(i),{16+i:1,i:-v})
assert len(cols)==29
S=s.Matrix.hstack(*[col for name,col in cols])
R=s.Matrix.vstack(base, v*base)
assert R.shape==(32,3) and S.shape==(32,29)
assert all(s.expand(z)==0 for z in R.T*S)
assert all(s.Poly(S[i,j],u,v).degree(u)<=1 and s.Poly(S[i,j],u,v).degree(v)<=1
           for i in range(32) for j in range(29))
print('PASS: 29 exact 32-row polynomial syzygies with degree <= (1,1)')
# The missing v^2 syzygy: an all-source Bob-infinity transfer, with no C^-1.
high=s.zeros(32,1)
for idx,coef in {27:v*v,20:-u*v,31:-u*v*v,16:-u*v,24:v,26:u*v}.items():
    high[idx]=coef
assert all(s.expand(z)==0 for z in R.T*high)
assert all(s.Poly(high[i],u,v).degree(u)<=1 and s.Poly(high[i],u,v).degree(v)<=2 for i in range(32))
print('PASS: exact additional u<=1,v<=2 syzygy transferring the full Bob-infinity charge')
uc,vc=s.symbols('u_c v_c',positive=True)
Rc=R.subs({u:uc,v:vc}); Sc=S.subs({u:uc,v:vc})
ref=[X0,Y0,C0]
assert s.factor(Rc[ref,:].det()) == uc
# All nonpivot original row coefficients of b are eliminated by the 13 first syzygies.
b=s.symbols('b0:32')
eta=[b[i]+vc*b[16+i] for i in range(16)]
sol=[]
for name,col in cols[:13]:
    ix={'Cv':6,'D':11,'m':4}.get(name)
    if ix is None: ix=int(name[2:])
    sol.append(eta[ix]/uc if ix==4 else (eta[ix]+vc*eta[4] if ix==9 else eta[ix]))
sol += list(b[16:])
res=s.Matrix(b)-Sc*s.Matrix(sol)
assert all(s.cancel(res[i])==0 for i in range(32) if i not in ref)
P=s.Matrix([res[i] for i in ref])
assert (s.simplify(Rc[ref,:].det()-uc)==0)
assert all(s.cancel(z)==0 for z in (Rc[ref,:].T*P-Rc.T*s.Matrix(b)))
assert all(s.simplify(z)==0 for z in Rc.T*Sc)
print('PASS: contact evaluation spans full 29-dimensional kernel, with explicit inverse using 1/u_c only')

# Source κ scalings and the two shifted (1+v²) copies.
Rhat=s.Matrix.vstack(s.diag(*kappa)*base,s.diag(*kappa)*(v*base))
q0=s.Matrix([kappa[i]*rows[i][1][3] for i in range(16)])
qhat=q0.col_join(v*q0)
Khat=s.diag(*(kappa+kappa))
assert all(s.simplify(z)==0 for z in Rhat.T*(Khat.inv()*S))
xi=s.Matrix([W,X,Y])
Fshort=sum((Rhat[i,0]*W+Rhat[i,1]*X+Rhat[i,2]*Y+qhat[i]*Z)**2 for i in range(32)) - 2*Du*Dv**2*Z**2
assert s.expand(Fshort-Dv*(F_expected-2*Du*Dv*Z**2))==0
print('PASS: exact (1+v²)*whole-parent 32-square-minus-one-square identity')

# A generic moving coefficient solution t = Khat^-1 S lambda, with 29 algebraic
# source numbers; all coefficients of every resulting state-linear root remain
# in the fixed (u<=1,v<=2) polynomial monomial subset of Z_(6,12).
qdeg=[s.Poly(s.expand(qhat[i]),u,v) for i in range(32)]
assert all(h.degree(u)<=1 and h.degree(v)<=2 for h in qdeg)
for i in range(32):
    for j in range(29):
        assert s.Poly(S[i,j],u,v).degree(u)<=1
        assert s.Poly(S[i,j],u,v).degree(v)<=1
print('PASS: source-scaled syzygy deformations preserve the K364 angular-degree bounds')
# For any original contact interpolation, subtract an interpolation of high(uc,vc)
# so the correction high-S*lambda_h vanishes at the contact. S*lambda_h
# has v-degree at most one, so it cannot affect the degree-four scalar pivot.
hi=s.Matrix([s.cancel(high[i].subs({u:uc,v:vc})) for i in range(32)])
hi_eta=[hi[i]+vc*hi[16+i] for i in range(16)]
hi_lam=[]
for name,col in cols[:13]:
    ix={'Cv':6,'D':11,'m':4}.get(name)
    if ix is None: ix=int(name[2:])
    hi_lam.append(hi_eta[ix]/uc if ix==4 else (hi_eta[ix]+vc*hi_eta[4] if ix==9 else hi_eta[ix]))
hi_lam += list(hi[16:])
assert all(s.cancel(z)==0 for z in hi-Sc*s.Matrix(hi_lam))
assert all(s.Poly(z,u,v).degree(v)<=1 for z in S*s.Matrix(hi_lam))
print('PASS: the high-degree transfer is compatible with every moving contact using source-only interpolation')
# Compute [v^4] of the scalar remainder for transfer tau*high.
# The general contact-interpolating component has degree <=1 in v.
# Therefore only the v^2 coefficients of high and q contribute.
tau=s.symbols('tau',real=True)
K2=[s.expand(k*k) for k in kappa]*2
q2=s.Matrix([s.Poly(qhat[i],v).coeff_monomial(v**2) for i in range(32)])
hi2=s.Matrix([s.Poly(high[i],v).coeff_monomial(v**2) for i in range(32)])
t2=s.Matrix([tau*hi2[i]/kappa[i%16] for i in range(32)])
v4=s.factor(-2*Du-2*(q2.T*t2)[0]-(t2.T*t2)[0])
want=s.factor(-2*(1+u*u)+2*tau*(1+u*u)-tau*tau*(1/(4*D+2)+u*u/2))
assert s.simplify(v4-want)==0
assert s.simplify(v4.subs(tau,2)-4*D/(2*D+1))==0
print('PASS: source-uniform positive [v^4] scalar pivot = 4D/(2D+1) when tau=2')

# New exact positive polynomial H in det B = 8 (1+u²)(1+v²) H(u²,v²).
BB=Rhat[:16,:].T*Rhat[:16,:]
H=s.cancel(BB.det()/(8*Du*Dv)); assert s.denom(H)==1
HP=s.Poly(s.expand(H),u,v)
expkeys={(4,4),(4,2),(4,0),(2,4),(2,2),(2,0),(0,4),(0,2)}
assert set(dict(HP.terms()))==expkeys
allmon=s.Poly(H,n,m,C,A,D,al,u,v)
for mon,coef in HP.terms():
    for _,scalar in s.Poly(s.expand(coef),n,m,C,A,D,al).terms():
        assert scalar>0,(mon,scalar)
print('PASS: det B/(8 Du Dv) has exactly eight strictly positive native-event coefficient polynomials')

output={
    'statement':'FULL-PARENT SYZYGY AND DET-B AUDIT; UNIVERSAL K364 ADMISSION UNPROVED',
    'base_row_order':[name for name,_,_ in rows],
    'base_row_factors':[str(k) for k in kappa],
    'base_row_shapes':[[str(a) for a in row] for _,row,_ in rows],
    'syzygy_order':[name for name,_ in cols],
    'extra_top_degree_syzygy':{str(i):str(high[i]) for i in range(32) if high[i]!=0},
    'top_bob_coefficient_tau':'-2*(1+u**2)+2*tau*(1+u**2)-tau**2*(1/(4*D+2)+u**2/2)',
    'top_bob_pivot_at_tau_2':'4*D/(2*D+1)>0 for every strict C2 source',
    'syzygy_sparse':[
       {'name':name,'entries':{str(i):str(col[i]) for i in range(32) if col[i]!=0}}
       for name,col in cols],
    'H_coefficients':{f'u^{ii}v^{jj}':str(s.factor(coef)) for ((ii,jj),coef) in HP.terms()},
    'interpolation': 'eta_i=b_i+v_c*b_(16+i), lambda_i=eta_i for nonpivot i except m; lambda_m=eta_m/u_c; lambda_shift_i=b_(16+i). Residual in pivot rows equals (R_pivot^T)^-1 R_c^T b.',
    'remaining': 'Find source-only syzygy gauges with scalar remainder -2DuDv²-2q^T t-t^T t an SOS in fixed Z_(6,12) Z features. The leading v^4 pivot is certified by high degree syzygy at tau=2, but the lower-degree scalar and contact-coupled PSD inequalities remain open. This is a sufficient ansatz, not an equivalent restriction of K364.'
}
Path(__file__).with_name('K364_SYZYGY_DATA.json').write_text(json.dumps(output,indent=2)+'\n')
print('ALL CHECKS PASSED; K364 strict C2 universal PSD admission is NOT proved')
