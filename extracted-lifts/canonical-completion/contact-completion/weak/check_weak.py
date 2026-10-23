"""Bounded independent weak recovery audit, with termwise rational clearing.

The received generator was read before execution; it only constructs symbolic
expressions and prints its44 records. No networking, file mutation or optimizer.
This checker does not certify source inequalities or Gram positive definiteness.
"""
from pathlib import Path
import sys
import contextlib, io, json, runpy, time
import sympy as s
HERE=Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()) as captured:
    g=runpy.run_path(str(HERE/'received_generator.py'))
(HERE/'received_records.txt').write_text(captured.getvalue(),encoding='utf8')
for name in ['B','T','j','p','H','n','mu','m','C','A','D','alpha','fc','qc','rc','gc','kc','dc','u','v','t','r','w','zeta','x','y','qz','rec','src']:
    globals()[name]=g[name]
half=s.Rational(1,2);h0=mu+half;kap=-kc
cb=n*j/(n+h0);Me=n*h0*j*j/(n+h0)
N=t-p*v*w;MM=-r-H*v*zeta;bc=v*w-u*zeta/j
q=v*zeta-T*u*w;R=zeta-u*v*w/B;Ru=v*(u*u-x)*w
Q2=q*q/(T*T*x);Rh2=R*R/qz;Bh2=bc*bc/y
Psi=(u*u/x+v*v/y)*w*w+(1+u*u*v*v/(x*y))*zeta*zeta/qz-4*u*v*w*zeta/(B*qz)
L=w-u*v*zeta/(B*qz)
E18=Q2+Rh2+v*v/y*L*L
W19=4*Psi+4*u*u/x*(Q2+Rh2)-v*v/y*(u*u/x-1)**2*w*w
Sigma=Me+fc/y;gamma=p*p/Sigma
tests=(u,v,t,r,w,zeta)
results={}
P=g['P'];E=g['E'];K=g['K'];Lsrc=g['L'];NN=g['NN'];FF=g['FF'];KK=g['KK']
# Explicitly cancel the common stationary denominator P BEFORE substitution.
# Every comparison below is affine in the listed source coefficients.
psrc={mu:j*p*(1+y)*(1+x)*(1+H*y)/2,
 n:NN/(2*(j-p)),C:Lsrc/(2*H*T),
 A:T*(1+x)*(1+H*y)*(K+H*j*y)/(2*(1-H*T)),
 D:(1+H*y)*(K+H*j*y)*(1-y*((1-H)*j-p))/(2*T*(1-H)),
 alpha:B**2*(K+H*j*y)*Lsrc/(2*E),
 fc:T*x*x*FF/2,
 rc:B*B*j*T*(1-p)*(1+x)*(1+H*y)*Lsrc/(2*E),
 gc:H*T*T*j*j*B**3*Lsrc/2,
 kc:-T*T*B*B*j*(1+x)*KK/(2*E)}
psrc[m]=psrc[mu]/(1-p)
psrc[qc]=psrc[A]*x-psrc[fc]
psrc[dc]=H*H*(psrc[A]+psrc[C]+psrc[alpha]+P/2)*y*y*qz-psrc[kc]
source_vars=list(psrc)
def native_times_P(c):
    out=0
    for powers,value in s.Poly(c,*source_vars).terms():
        assert sum(powers)<=1,(powers,c)
        if sum(powers)==0:out+=value*P
        else:out+=value*psrc[source_vars[powers.index(1)]]
    return out

def zero_coeffs(expr,name,substitution=None):
    """Clear each testing coefficient separately; never expand the native parent."""
    poly=s.Poly(s.expand(expr),*tests)
    failures=[];nontrivial=0
    for powers,c in poly.terms():
        if c==0:continue
        nontrivial+=1
        if substitution is not None:c=native_times_P(c)
        num,den=s.fraction(s.cancel(c))
        if num!=0:
            failures.append((list(powers),str(s.factor(num))))
        if substitution is not None and nontrivial%10==0:
            print('  ',name,nontrivial,'cleared',flush=True)
    results[name]={'testing_coefficients_cleared':nontrivial,'zero':not failures,'failures':failures}
    print(name,nontrivial,'coefficients',len(failures),'failures',flush=True)
    return not failures

# Collapse only the identical6/7/8 diagonal records, using their exact circuit
# sum of squares. Every other actual record is contracted from the given list.
normal_f={i+1:g['f'][i] for i in range(11) if i not in [5,6,7]}
contract=0
for ii,jj,eu,ev,c in rec:
    if ii in [6,7,8]:
        assert ii==jj
        if ii==6:
            assert [z for z in rec if z[:4]==(7,7,eu,ev)][0][4]==c
            assert [z for z in rec if z[:4]==(8,8,eu,ev)][0][4]==c
            contract+=u**eu*v**ev*c*fc*Psi
        continue
    term=c*normal_f[ii]*normal_f[jj]
    # All square roots involved here have nonnegative pure-source arguments.
    term=s.powdenest(s.powsimp(s.expand_power_base(term,force=True),force=True),force=True)
    term=s.simplify(term)
    assert not any(power.exp.q!=1 for power in term.atoms(s.Pow)),(ii,jj,term)
    contract+=(1 if ii==jj else 2)*u**eu*v**ev*term

RP=(g['mu']*p/(1-p)*N*N+(mu*p/(1-p)+h0+alpha)*u*u*N*N+half*v*v*N*N
    +n*v*v*(N+p*bc)**2+(C+D+half)*u*u*MM*MM
    +(C+A+alpha+half)*v*v*MM*MM+C*u*u*v*v*MM*MM
    +dc*v*v/y*Q2+half*p*p*v*v*bc*bc+gc*v*v/y*L*L
    +kap*(R+v*q/(T*B))**2/qz)
certificate=contract+RP+16*x*x*y*W19+8*gamma*fc*y*v*v*E18+gamma*T*fc*y*u*u*E18
Fs=(qc*Q2+n*j*j*bc*bc+rc*Rh2+dc*v*v/y*Q2
    +(n*p*p*y*y+half*p*p*y*y-gc)*v*v/y*Bh2
    +fc*Psi+gc*v*v/y*Psi+kc*zeta*zeta/qz*(v*v/y-1)**2)
Kc=g['Kc'];KR=g['KR']
LN=n*j*bc+n*p*v*v*bc-j*alpha*bc+half*p*v*v*bc+KR*Ru+(Kc/j)*u*v*q
Cq=A/T+v*v*(C*H*(1+u*u)+A*H+(alpha+half)*H)
LM=Cq*q+mu*u*v*bc+v*(alpha-B*T*C*H*u*u)*R
HN=h0*(1+u*u)+half*v*v+alpha*u*u+n*(1+v*v)
HM=(1+u*u)*(C*(1+v*v)+D+half)+A*(1+v*v)+(alpha+half)*v*v
parent=Fs+2*LN*N+2*LM*MM+HN*N*N+HM*MM*MM+2*(alpha+half)*u*v*N*MM+mu*p/(1-p)*(1+u*u)*N*N
zero_coeffs(certificate-parent,'all44_record_restoration_to_RFN_parent')
zero_coeffs(gamma*Me*v**4*bc**2+gamma*fc*(y*v*v+v**4)*Psi-p*p*v**4*bc**2
            -gamma*fc*y*v*v*E18,'purchase1')
zero_coeffs(gamma*T*Me*u*u*v*v*bc**2+gamma*T*fc*u*u*(y+v*v)*Psi-T*p*p*u*u*v*v*bc**2
            -gamma*T*fc*y*u*u*E18,'purchase2')

# Original complete quarter-parent, with independent normal coordinates.
W=w-u*r;X=t-u*zeta;Y=r+u*w;Z=zeta+u*t;Du=1+u*u;Dv=1+v*v
original=(m*(v*W-X+u*v*Y-u*Z)**2+n*Du*Dv*X*X+C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y
          +(D+half)*(u*W+u*v*X-Y-v*Z)**2+(alpha+half)*Du*(v*Y-Z)**2
          +half*Dv*(X+u*Z)**2-half*Du*Dv*Z*Z)
# Original=Du*parent. Source substitution is applied only after extracting
# each small testing coefficient, so denominators never grow across state rows.
zero_coeffs(original-Du*parent,'full_original_parent_native_dictionary',src)

# Independently compare the scalar dictionary to its coefficient definitions.
K0=C*H*H+(D+half)*(1-H)**2
Ma=A*x;Mb=(mu*(1-p)+(n+half)*p*p)*y;Mc=alpha*qz
Mab=(A+alpha*(1-p)**2+mu*(1-p)+half*(1-p)**2)*x*y
Mac=n*x*qz;Mbc=(K0+A*H*H+2*H*alpha+H-half)*y*qz
Mabc=(K0+n)*x*y*qz;Mab2c=C*H*H*x*y*y*qz
definitions={kc:(Ma-Mb+Mac+Mab2c-Mbc)/2,
 rc:Mc-Mb+Mac-kc,fc:Mb-Mac,qc:Ma-fc,
 dc:H*H*(A+C+alpha+half)*y*y*qz-kc,gc:Mab2c}
for symbol,value in definitions.items():
    zero_coeffs(symbol-value,'source_mass_'+str(symbol),src)

# Direct original GE1 -> compact dictionary, with the known common denominator
# e=1+Bp removed first. These polynomial identities establish the six actual
# event weights, not only the internally consistent scalar mass definitions.
e=1+B*p;boN=B*(1+H*T*x);hoN=B*(1-H*T);rN=j-p
deltaN=(x*e*e+(1+x)*boN*e+(1+x)*boN*boN
        +2*x*(1+x)*boN*rN/j+x**3*rN*rN/(j*j)+2*x*x*rN*e/j
        +y*((1+x)*boN*e+x*e*e+x*x*(2*rN*e/j-rN*rN/(j*j))))
MN=(1+x)*B*(K+H*j*y)*(hoN+y*e)-x*(1+y)*p*H*y*(1+x)**2/j
GE1identities={
 'GE1_delta_times_e_squared':deltaN-B*(1+x)*P,
 'GE1_M_times_e_squared':MN-B*B*(1+x)*Lsrc,
 'GE1_G_times_e':boN+x*(hoN+rN/j-e)-B*(1-p)*(1+x),
 'GE1_n_numerator':j*(1+x)*(K+H*j*y)*(1+H*y)-(j-p)*P-NN,
 'GE1_B_plus_zr_times_e':boN+x*rN/j-B*(K+H*j*y),
 'GE1_B_minus_zyr_times_e':boN-x*y*rN/j-B*(1-y*((1-H)*j-p)),
 'GE1_one_minus_Ho_times_e':e-hoN-E}
for name,expr in GE1identities.items():zero_coeffs(expr,name)

# Literal original rational source substitution (no implicit contact equation).
z0=x;y0=y;r0=rN/(j*e);H0=hoN/e;B0=boN/e
cs0=(1+x)/e;t0=H*cs0;L0=y0*t0;gap0=(1-H)*cs0
delta0=B*(1+x)*P/e**2;M0src=B**2*(1+x)*Lsrc/e**2
G0=B*(1-p)*(1+x)/e;S0=cs0*(1+H*y)
native_GE1={n:((B0+x*r0)*cs0*S0-r0*delta0)/(2*r0*delta0),
 m:x*x*(1+y)*(1-r0)*cs0*S0/(2*G0*delta0),
 C:cs0*M0src/(2*L0*delta0),
 A:y*cs0*(B0+x*r0)*S0/(2*H0*delta0),
 D:(1+x)*(B0+x*r0)*(B0-x*y*r0)*S0/(2*y*gap0*delta0),
 alpha:(B0+x*r0)*M0src/(2*(1-H0)*delta0)}
for symbol,value in native_GE1.items():
    num,_=s.fraction(s.cancel(P*value-psrc[symbol]))
    results['literal_GE1_'+str(symbol)]={'zero':num==0}
    assert num==0,symbol

(HERE/'records.json').write_text(json.dumps([[ii,jj,eu,ev,str(c)] for ii,jj,eu,ev,c in rec],indent=2)+'\n',encoding='utf8')

(HERE/'identity_checks.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8')
assert all(r['zero'] for r in results.values()),'See identity_checks.json for exact discrepancies.'
print('All requested rational identities passed.',flush=True)
