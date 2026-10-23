"""Exact independent reconstruction of the full-R1 thirteen-row normal core.

The completed Pro reply uses the whole Jalpha row, avoiding duplicate gA/Ru
payments. This script checks its full RFN-parent interface, covariance,
two purchase identities, Schur arithmetic and explicit source-box margins.
No optimizer, sampling or received pickle is used. This is not a compiler
or an admission proof for the complementary testing regions.
"""
from pathlib import Path
import contextlib,io,json,runpy
from fractions import Fraction as F
import sympy as s
HERE=Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    g=runpy.run_path(str(HERE/'source_dictionary.py'))
for name in ['B','T','j','p','H','n','mu','m','C','A','D','alpha','fc','qc','rc','gc','kc','dc','u','v','t','r','w','zeta','x','y','qz']:
    globals()[name]=g[name]
half=s.Rational(1,2);h0=mu+half;kap=-kc
Me=n*h0*j*j/(n+h0);Sigma=Me+fc/y;theta=Me/Sigma
phi=s.symbols('phi',real=True)
aN,bN,tau,rhoQ,KR,Kc,aM,hH=s.symbols('aN bN tau rhoQ KR Kc aM hH',real=True)
N=t-p*v*w;MM=-r-H*v*zeta;bc=v*w-u*zeta/j
q=v*zeta-T*u*w;Rphys=zeta-u*v*w/B;Ru=v*(u*u-x)*w
Jalpha=j*alpha*bc;L=alpha-B*T*C*H*u*u
O=s.sqrt(T)*v*L*Rphys;H1=p*v*v*bc;H2=s.sqrt(T)*p*u*v*bc
Q2=q*q/(T*T*x);R2=Rphys*Rphys/qz;Bs2=bc*bc/y
Psi=(u*u/x+v*v/y)*w*w+(1+u*u*v*v/(x*y))*zeta*zeta/qz-4*u*v*w*zeta/(B*qz)
closed=w-u*v*zeta/(B*qz)
tests=(u,v,t,r,w,zeta)
checks={}
def zero(expr,name):
    failures=[];count=0
    for powers,c in s.Poly(s.expand(expr),*tests).terms():
        if c==0:continue
        count+=1
        num,_=s.fraction(s.cancel(c))
        if num!=0:failures.append([list(powers),str(s.factor(num))])
    checks[name]={'zero':not failures,'coefficients_cleared':count,'failures':failures}
    assert not failures,(name,failures[:3])
    print(name,'passed',flush=True)

# Actual symmetric core entries, zero-based row order given in the result.
K=s.diag(1,half,16,16,8,aM,16,8,half,half,s.Rational(1,4),s.Rational(1,4),s.Rational(1,4))
def entry(i,j,c):K[i,j]=K[j,i]=c
for target,coefficient in [(2,KR),(3,-1),(4,half)]:
    entry(0,target,aN*coefficient);entry(1,target,-bN*coefficient)
entry(5,6,hH);entry(5,7,hH*mu/p)
entry(0,8,aN*(Kc/j)*rhoQ*u*v);entry(1,8,-bN*(Kc/j)*rhoQ*u*v)
entry(0,5,aN*(alpha+half)*tau*u*v);entry(1,5,-bN*(alpha+half)*tau*u*v)
Cq=A/T+H*v*v*(C*(1+u*u)+A+alpha+half)
entry(5,8,tau*rhoQ*Cq)
records=[]
for i in range(13):
    for k in range(i,13):
        for powers,c in s.Poly(K[i,k],u,v).terms():
            if c!=0:records.append([i,k,list(powers),str(s.factor(c))])
assert len(records)==28
assert sorted({tuple(r[2]) for r in records})==[(0,0),(0,2),(1,1),(2,2)]

# Exact polynomial core decomposition against the complete RFN parent.
# The current identity uses N0=aN*Nhat-bN*E and
# Nhat^2+E^2=h0*N0^2+n*(N0+j*b)^2, with E^2=Me*b^2.
KRactual=g['KR'];Kcactual=g['Kc']
Tcost=16*Ru*Ru+16*Jalpha*Jalpha+16*T*v*v*L*L*Rphys*Rphys+8*H1*H1+8*T*p*p*u*u*v*v*bc*bc
core=(h0*N*N+n*(N+j*bc)**2-half*Me*bc*bc+(C+D+half+A)*MM*MM
      +half*qc*Q2+half*rc*R2+s.Rational(1,4)*fc*Psi+Tcost
      +2*N*(KRactual*Ru-Jalpha+H1/2)
      +2*v*L*MM*Rphys+2*mu*u*v*MM*bc
      +2*u*v*N*((Kcactual/j)*q+(alpha+half)*MM)+2*Cq*MM*q)
mp=mu*p/(1-p)
positive=(mp*(1+u*u)*N*N+(h0+alpha)*u*u*N*N+half*v*v*N*N
          +n*v*v*(N+p*bc)**2+(C+D+half)*u*u*MM*MM
          +(C+A+alpha+half)*v*v*MM*MM+C*u*u*v*v*MM*MM
          +dc*v*v/y*Q2+half*p*p*v*v*bc*bc+gc*v*v/y*closed*closed
          +kap*(Rphys+v*q/(T*B))**2/qz+gc*(Q2+R2))
Fs=(qc*Q2+n*j*j*bc*bc+rc*R2+dc*v*v/y*Q2
    +(n*p*p*y*y+half*p*p*y*y-gc)*v*v/y*Bs2
    +fc*Psi+gc*v*v/y*Psi+kc*zeta*zeta/qz*(v*v/y-1)**2)
LN=n*j*bc+n*p*v*v*bc-j*alpha*bc+half*p*v*v*bc+KRactual*Ru+(Kcactual/j)*u*v*q
LM=Cq*q+mu*u*v*bc+v*L*Rphys
HN=h0*(1+u*u)+half*v*v+alpha*u*u+n*(1+v*v)
HM=(1+u*u)*(C*(1+v*v)+D+half)+A*(1+v*v)+(alpha+half)*v*v
parent=(Fs+2*LN*N+2*LM*MM+HN*N*N+HM*MM*MM+2*(alpha+half)*u*v*N*MM
        +mp*(1+u*u)*N*N)
residual=(half*Me*bc*bc+half*qc*Q2+half*rc*R2+s.Rational(3,4)*fc*Psi
          -gc*Psi-2*kap*v*v/y*Q2-2*kap*R2-Tcost)
zero(Jalpha+alpha*Ru/B+alpha*u*Rphys,'whole_Jalpha_identity')
zero(parent-core-positive-residual,'entire_RFN_parent_equals_core_plus_retained_plus_scalar_residual')

# Both RFN18 forms are valid: their common identity is
# (1+U^2)Psi=Rstar^2+Bstar^2+U^2 Qstar^2+U^2 closed^2.
A27=(1+u*u/x)*Psi-R2
B29=2*(Psi+v*v/y*Q2+u*u/x*Bs2)-R2
zero(A27-Bs2-u*u/x*Q2-u*u/x*closed*closed,'RFN18_original_Rstar_subtraction')
zero((1+u*u/x)*Psi-Bs2-R2-u*u/x*Q2-u*u/x*closed*closed,'RFN18_Pro_Bstar_subtraction')
A1=alpha*Ru/B;B1=alpha*u*Rphys
purchaseJ=(2*(1-phi)*A1*A1+2*(1-phi)*B1*B1+phi*j*j*alpha*alpha*y*(1+u*u/x)*Psi-Jalpha*Jalpha)
remainderJ=(1-phi)*(A1-B1)**2+phi*j*j*alpha*alpha*y*(R2+u*u/x*Q2+u*u/x*closed*closed)
zero(purchaseJ-remainderJ,'Pro_equation1_exact_whole_Jalpha_purchase')
W19=4*Psi+4*u*u/x*(Q2+R2)-v*v/y*(u*u/x-1)**2*w*w
RA=W19+4*u*u/x*A27
RB=8*u*u/x*Psi+8*u**4/x**2*Q2+(8*u**4/x**2+2*(u*u/x-1)**2)*Bs2-v*v/y*(u*u/x-1)**2*w*w
s0=u*u
rPsi=(1-phi)*4*x*x*y/fc+phi*(8*theta*x*y*s0/fc+4*(x*x+x*s0+s0*s0)/Sigma)
rQ=(1-phi)*4*x*y*s0/qc+phi*(8*theta*y*s0*s0+4*(1-theta)*x*y*s0)/qc
rE=phi*(8*s0*s0+2*(s0-x)**2)/Sigma
rR=(1-phi)*4*x*y*s0/rc
zero(fc*rPsi*Psi+qc*rQ*Q2+Me*rE*bc*bc+rc*rR*R2-Ru*Ru
     -x*x*y*((1-phi)*W19+phi*(1-theta)*RA+phi*theta*RB),'Pro_equation2_exact_smooth_radial_purchase')

# Covariance of all nine original source columns after the state rotation.
W=w-u*r; X=t-u*zeta; Y=r+u*w; Z=zeta+u*t; Du=1+u*u; Dv=1+v*v
columns = {'lambda':(W+v*X)**2, 'n':Dv*X*X, 'm':Du*(v*w-t)**2,
           'C':Du*Dv*r*r, 'A':Dv*Y*Y, 'D':Du*(r+v*zeta)**2,
           'alpha':(v*Y-Z)**2, 'beta':Du*Dv*t*t}
columns['d']=(columns['D']+columns['alpha']+columns['beta']-Dv*Z*Z)/2
cov={}
for a0 in [-1,1]:
    for b0 in [-1,1]:
        signs=[a0,b0,1,a0*b0,b0,a0]
        cov[str((a0,b0))]=all(all(c==c*__import__('math').prod(si**ei for si,ei in zip(signs,mon))
            for mon,c in s.Poly(poly,*tests).terms()) for poly in columns.values())
assert all(cov.values())
checks['all_nine_source_covariance']=cov

# Independent Schur arithmetic and uniform perturbation estimates.
schurN=2*(F(32,255)+F(4,127));schurM=F(4,255)+F(16,1143)
assert schurN==F(10168,32385)<F(7,16)
assert schurM<F(1,16)
delta=F(1,1024);ell=1+delta+delta**2+delta**3
f0=1-delta-delta**3;Rdelta=1+2*delta+2*delta**2+4*delta**3+delta**4*(1+2*delta)
assert (1+delta**2)*ell/(2*(1-delta))<1 # alpha<=B(B+x)
assert (1+delta)*(1+delta**2)**2/2<F(2,3) # mu/p
assert (1+delta)*(1+delta**2)**2/(2*(1-delta))<F(2,3) # A/T
assert 1+4*delta+2*delta**2+6*delta**3+delta**4*(1+delta)<2 # NN/j
assert (1-2*delta)/(4*Rdelta)>F(1,8) # qc/(Tx)
assert f0/(4*Rdelta)>F(1,8)
aMupper=max(ell,(1+delta**2)**3)/2+delta/8+delta**2/6
assert aMupper<F(3,4)
assert F(2,3)+2*delta**2+2*delta<1
# Operator norm bounds: K11^2<=578T/9; K02,K22 use squared rational bounds.
assert F(578,9)*delta<1
K02square=2*(ell/2+delta*(2*delta/3+2*delta**2+F(1,2)))**2
assert K02square<1
assert 2*(ell/2)**2<1
assert F(8,9)<1 # constant perturbation<=2sqrt(2)T/3<T
gatebudget=F(1,2**37)+F(1,2**50)+F(1,2**74)
assert F(1,16)-delta-gatebudget>F(1,32)
checks['Schur_and_margin_exact']={'N_correction':str(schurN),'M_correction':str(schurM),
 'aM_lower':str(f0/(4*Rdelta)),'aM_upper':str(aMupper),
 'source_KR_upper':str(F(2,3)+2*delta**2+2*delta),
 'K11_norm_squared_upper':str(F(578,9)*delta),'K02_norm_squared_upper':str(K02square),
 'reserve_lower':str(F(1,16)-delta-gatebudget)}
out={'row_order':['Nhat','Ecurrent','Ru','Jalpha','H1','Mhat','O','H2','Qhat','Rhat','psi1','psi2','psi3'],
 'indexing':'zero-based actual symmetric entries; off-diagonal factor2 in contraction',
 'records':records,'source_definitions':{'aN':'1/sqrt(n+mu+1/2)','bN':'sqrt(n/((mu+1/2)*(n+mu+1/2)))',
 'tau':'sqrt(T*H*(1-H))','hH':'sqrt(H*(1-H))','rhoQ':'T*sqrt(x/qc)',
 'KR':str(g['KR']),'Kc':str(g['Kc']),'aM':'T*H*(1-H)*(C+D+1/2+A)'},
 'scope':'normal core and exact scalar residual only; simultaneous purchases separately audited'}
(HERE/'normal_core_records.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
(HERE/'normal_core_checks.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf8',newline='\n')
print('28 exact core records, full-parent interface, covariance and Schur margins passed.')
