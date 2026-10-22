"""Portable exact C1 feature extraction. Requires Python 3 and SymPy only.

Checks identities introduced by the fixed-feature compilation; it does not
claim a global source/chart cover or replace the physical parent audit.
"""
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = {}
v, w = sp.symbols('v w')
u = v*v
names = ('nu R0 p D m q z sig k dG rho S q0 c a0 L0 b Q '
         'rdG rrho rsig rlambda rq0 rk KA a ell C '
         'eyu eyw er0 eru erw ec0 ecu ecw Rx Bx Ay chi').split()
src = dict(zip(names,sp.symbols(' '.join(names))))
globals().update(src)
N, JC, JL, GG, AT, TT, d1,d2,d3,d4,UU,KK = roots = sp.symbols(
    'Nx JC JL GammaStar Atheta Ttheta d1 d2 d3 d4 Ustar Kstar')
ny = Q*JC+b*JL
rr = JC-JL
f1 = rdG*u*GG/S+rrho*(k*d1-u*d2+(u-k)*d3)/S-rho*rsig*TT/S
f2 = (R0*rdG*q*w*GG/S+R0*rho*sig*z*rr/S
      -R0*rrho*(q*w+z*rsig)*(d2-d3)/S
      +(R0*q0*z/S+R0*rlambda*rrho*q*w/S)*TT-8*p*q*u*w*f1)
epsC = 2/sp.Symbol('lam')+u+k
epsG = 2*a0*u+(u+k)/L0
ecommon = ((eyu*u*u+eyw*w*w)*ny*ny
    +(er0+eru*u*u+erw*w*w)*rr*rr
    +(ec0+ecu*u*u+ecw*w*w)*ny*rr-4*nu*N*f1+JC*f2
    -4*u*q*w*(2*m+2*D+1)*N*ny+u*q*w*Rx*N*rr
    +q*w*Bx*N*AT+u*Ay*ny*AT+chi*u*AT*rr)
weighted = a0*(u*UU*UU+2*rk*v*UU*KK+k*KK*KK)
ecorr = (ecommon-rho*sig*rr*rr/2+rq0*v*UU*rr
    -c*dG*epsC*(d1*d1+d2*d2+d3*d3+d4*d4)/S
    -c*dG*u*u*rr*rr/S-c*epsG*GG*GG/S-c*weighted/S)
diag = [2,1,1,1,sp.Rational(1,2),sp.Rational(1,2),
        sp.Rational(1,2),sp.Rational(1,2),sp.Rational(1,2),
        sp.Rational(1,2),1,sp.Rational(1,2)]
form = sp.expand(sum(a*r*r for a,r in zip(diag,roots))+ecorr)
gram = sp.hessian(form,roots)/2
support = sorted({mon for i in range(12) for j in range(i,12)
                  for mon,_ in sp.Poly(gram[i,j],v,w).terms()})
nonconstant = [(0,1),(0,2),(1,0),(2,0),(2,1),(4,0),(4,1)]
assert support == [(0,0)]+nonconstant
checks['main_testing_support_exact'] = support

# Print every fixed direction and its source-constant scalar coefficients.
directions=[]
channels=[]
rebuilt=sp.zeros(12)
def channel(vec,scalar,kind):
    ident=f'D{len(directions)+1:03d}'
    directions.append({'id':ident,'kind':kind,'vector':vec})
    poly=sp.Poly(sp.expand(scalar),v,w)
    channels.append({'direction':ident,
        'source_constant':str(poly.coeff_monomial(1)),
        'source_coefficients':[str(poly.coeff_monomial(v**a*w**bb))
                               for a,bb in nonconstant]})
    for i,vi in enumerate(vec):
        if vi:
            for j,vj in enumerate(vec):
                if vj: rebuilt[i,j] += vi*vj*scalar
for i in range(12):
    vec=[0]*12;vec[i]=1
    channel(vec,gram[i,i]-sp.Rational(11,48),'coordinate')
for i in range(12):
    for j in range(i+1,12):
        for sign in (1,-1):
            vec=[0]*12;vec[i]=1;vec[j]=sign
            channel(vec,(sp.Rational(1,48)+sign*gram[i,j])/2,
                    'pair_plus' if sign==1 else 'pair_minus')
assert len(directions)==144
assert all(sp.expand(rebuilt[i,j]-gram[i,j])==0
           for i in range(12) for j in range(12))
checks['main_144_direction_identity_exact']=True

# Normal identities, using positive source square roots l,s,k,q,d.
l,s,kc,qc,dc,uc,wc,XX,YY,ZZ,WW,gcommon,Rc,pc = sp.symbols(
    'l s kc qc dc uc wc XX YY ZZ WW gcommon Rc pc')
rho_c=l*l*kc*kc/(s*s)
Sc=dc*dc+rho_c
B=[l*s*uc*ZZ,l*kc*uc*wc*WW/s,l*kc*kc*qc*WW,l*kc*wc*ZZ/qc]
dd=[bb-gcommon for bb in B]
gamma=s*s*ZZ-uc*wc*WW
rrc=wc*ZZ/qc-uc*qc*WW
tc=kc*qc*(uc-kc)*WW
f1c=(dc*uc*(dc*gamma)/Sc+(l*kc/s)*(kc*dd[0]-uc*dd[1]
    +(uc-kc)*dd[2])/Sc-rho_c*s*tc/Sc)
assert sp.cancel(f1c-uc*gamma)==0
f2c=(Rc*dc*qc*wc*(dc*gamma)/Sc+Rc*rho_c*s*s*qc*qc*rrc/Sc
    -Rc*(l*kc/s)*(qc*wc+qc*qc*s)*(dd[1]-dd[2])/Sc
    +(Rc*l*l*kc*qc*qc/Sc+Rc*l*(l*kc/s)*qc*wc/Sc)*tc
    -8*pc*qc*uc*wc*f1c)
assert sp.cancel(f2c-qc*wc*(Rc-8*pc*uc*uc)*gamma)==0
checks['both_full_normal_identities_exact']=True

# Averaged scalar identity and weighted angular identity.
FA=uc*(ZZ-kc*wc*WW/(s*s))
AA=qc*(uc*uc-kc*kc)*WW
Pbeta=uc*ZZ-kc*kc*wc*WW/(s*s)
Sorig=l*l*s*s*Pbeta**2+rho_c*(uc*uc-kc*kc)*(wc*wc-s*s*qc*qc)*WW**2
U2=l*l*kc*uc*qc*qc*(uc-kc)**2*WW**2
K2=l*l*kc*kc*qc*qc*(uc-kc)**2*WW**2
cross=l*l*kc*uc*qc*(uc-kc)*WW*rrc
assert sp.cancel(Sorig-(l*l*s*s*FA**2+2*U2+K2+2*cross))==0
weighted_c=(kc/(s*s))*l*l*kc*qc*qc*(uc-kc)**2*WW**2*(uc+kc)**2
assert sp.cancel(weighted_c-rho_c*AA**2)==0
checks['averaged_scalar_and_weighted_angular_identities_exact']=True

# Whole-parent assembly: its sole formal residual is the verified weighted
# angular relation, before evaluating the actual dependent state roots.
CP,FAformal=sp.symbols('Cplus FA')
cc=sum(d*d for d in (d1,d2,d3,d4))+CP
old=(ecommon+GG*GG+(2+(2+4*m)*w*w)*N*N
    +(4*a+2*z*(2*C+1))*JC*JC+4*ell*JL*JL+KA*AT*AT
    +L0*FAformal**2/2+UU*UU+KK*KK/2+rq0*v*UU*rr+cc/2
    -rho*sig*rr*rr/2)
ordinary=((2+4*m)*w*w*N*N+(4*a+2*z*(2*C+1)-1)*JC*JC
    +(4*ell-1)*JL*JL+L0*FAformal**2/2
    +(KA+c-1)*(AT*AT+TT*TT)/2+(KA+c)*(AT*AT-TT*TT)/2)
budget=(c*dG*epsC*cc/S+c*dG*u*u*rr*rr/S+c*epsG*GG*GG/S
        -c*dG*AT*AT/S)
qC=sp.Rational(1,2)-c*dG*epsC/S
assembly=sp.expand(form+ordinary+budget+qC*CP-old
                  -c*(rho*AT*AT-weighted)/S)
assert sp.cancel(assembly.subs(S,dG+rho))==0
checks['whole_parent_assembly_mod_verified_relation_exact']=True

# Eight angular squares with every state denominator already cleared.
au0,av0,az0,aa0,mm,du,dw,dW,dZ = sp.symbols(
    'au0 av0 az0 aa0 mm du dw dW dZ')
au,av,az,astar=au0**2,av0**2,az0**2,aa0**2
ut,wt,Wt,Zt=du**2,dw**2,dW**2,dZ**2
rad=dw*dW*dZ
af=[mm*(au*av*az*ut*wt*Zt-Wt),
    sp.sqrt(2)*mm*(av*az*wt*Zt-au*ut*Wt),
    mm*((aa0+1/aa0)*au*av*ut*wt*Wt-aa0*au*az*ut*Zt-az*Zt/aa0)]
af.append(au0*du*af[2])
af.append(aa0*mm*(Wt-av*az*wt*Zt))
af.append(au0*du*af[4])
af.append(sp.sqrt(2)*aa0*av0*az0*mm*(au*ut-1)*rad)
af.append(au0*du*af[6])
def angular_lhs(sign):
    zt=sign*Zt
    psi=mm**2*(Wt**2+au**2*av**2*ut**2*wt**2*Wt**2
        +au**2*az**2*ut**2*zt**2+av**2*az**2*wt**2*zt**2
        -4*au*av*az*ut*wt*Wt*zt)
    mq=mm*(av*az*wt*zt-au*ut*Wt)
    mr=mm*(az*zt-au*av*ut*wt*Wt)
    mt=mm*Wt*((au*ut)**2-1)
    return ((2+astar*(au*ut+1))*psi+(au*ut)**2*mq*mq
       +(2*au*ut+(au*ut+1)/astar)*mr*mr-mt*mt)
assert sp.expand((angular_lhs(1)-sum(ff*ff for ff in af))*aa0**2)==0
cx2=4;cx1=8*astar+8+4/astar;cx0=16+8*astar+4/astar
negative=mm**2*au*av*az*ut*wt*Wt*Zt*(cx2*(au*ut)**2+cx1*au*ut+cx0)
assert sp.expand((angular_lhs(-1)-angular_lhs(1)-negative)*aa0**2)==0
checks['eight_cleared_angular_roots_exact']=True
checks['three_negative_state_angular_roots_exact']=True

# Circuit identity and its six testing monomials, using fourth-power
# source/testing variables so every calculation is polynomial.
a1,a2,a3,a4,tu,tw,tW,tZ=sp.symbols('a1 a2 a3 a4 tu tw tW tZ')
cs=[a1**4,a2**4,a3**4,a4**4]
ct=[tu**3*tw*tW*tZ**3,tu**3*tw**3*tW**3*tZ,
    tu*tw*tW**3*tZ,tu*tw**3*tW*tZ**3,
    tu**3*tw**2*tW**2*tZ**2,tu*tw**2*tW**2*tZ**2]
cg=a1*a2*a3*a4
cm=[tu**4*tZ**4,tu**4*tw**4*tW**4,tW**4,tw**4*tZ**4]
cG=cg*tu**2*tw**2*tW**2*tZ**2
cp=2*cG*sum(si*mi for si,mi in zip(cs,cm))-8*cG*cG
cp_roots_squared=(2*cg*(a1*a1*ct[0]-a2*a2*ct[1])**2
    +2*cg*(a3*a3*ct[2]-a4*a4*ct[3])**2
    +4*cg*(a1*a2*ct[4]-a3*a4*ct[5])**2)
assert sp.expand(cp-cp_roots_squared)==0
delta,ac,qc0=sp.symbols('delta ac qc0')
assert sp.expand(qc0-ac*uc-((qc0-ac*delta)+ac*(delta-uc)))==0
checks['circuit_six_positive_features_identity_exact']=True
checks['circuit_negative_feature_and_scalar_split_exact']=True

# All 220 labelled testing features.
h=['X','Y','Z','W','u*Z','w*Z','u*W','u*w*W','u^2*W',
   'u^(3/2)*W','u^(1/2)*W','sqrt(u*w*W*Zplus)']
bt=['w','w^2','sqrt(u)','u','u*w','u^2','u^2*w']
factors=['1']+[f'sqrt(gamma7{sgn}({f}))' for f in bt for sgn in ('+','-')]
features=[]
for i,hh in enumerate(h,1):
    for j,fac in enumerate(factors):
        features.append({'id':f'M{i:02d}_{j:02d}','category':'main',
                         'expression':hh if j==0 else f'({hh})*({fac})'})
ordinary_features=['w*X','u*w*Z','u*X','u^2*w*W','u*Y','w*Y','u*w*Y']
for i,ff in enumerate(ordinary_features,1):
    features.append({'id':f'O{i:02d}','category':'ordinary','expression':ff})
angular_features=['W','u*w*abs(Z)','w*abs(Z)','u*W','u*w*W','u*abs(Z)',
    'abs(Z)','u^(3/2)*w*W','u^(3/2)*abs(Z)','sqrt(u)*abs(Z)',
    'sqrt(u)*W','sqrt(u)*w*abs(Z)','u*rstar','rstar','u^(3/2)*rstar',
    'sqrt(u)*rstar','Inegative*u^(3/2)*rstar','Inegative*u*rstar',
    'Inegative*sqrt(u)*rstar']
for i,ff in enumerate(angular_features,1):
    features.append({'id':f'A{i:02d}','category':'angular','expression':ff})
circuit_features=['u^(3/4)*w^(1/4)*W^(1/4)*Zplus^(3/4)',
    'u^(3/4)*w^(3/4)*W^(3/4)*Zplus^(1/4)',
    'u^(1/4)*w^(1/4)*W^(3/4)*Zplus^(1/4)',
    'u^(1/4)*w^(3/4)*W^(1/4)*Zplus^(3/4)',
    'u^(3/4)*w^(1/2)*W^(1/2)*Zplus^(1/2)',
    'u^(1/4)*w^(1/2)*W^(1/2)*Zplus^(1/2)','sqrt(u*w*W*Zminus)']
for j in (0,1):
    for i,ff in enumerate(circuit_features,1):
        features.append({'id':f'C{i:02d}_{j}','category':'circuit',
          'expression':ff if j==0 else f'({ff})*sqrt(delta-u)'})
assert len(features)==220 and len({f['id'] for f in features})==220
checks['feature_counts']={'main':180,'ordinary':7,'angular':19,'circuit':14,
                         'total':220}
checks['square_root_counts']={'main':144*15,'ordinary':17,'angular':11,
                             'circuit':8,'total':2196}

# Source-linear coordinate expressions of the twelve main roots.
main_coordinates=[{'M01_00':'1','M05_00':'h'},
 {'M02_00':'1/q','M07_00':'-b*q'},
 {'M02_00':'1/q','M06_00':'-1/q','M07_00':'Q*q'},
 {'M03_00':'rdG*sig','M08_00':'-rdG'},
 {'M09_00':'q','M04_00':'-q*k^2'},
 {'M07_00':'k*q','M04_00':'-k^2*q'},
 {'M05_00':'rlambda*rsig','M12_00':'-rlambda*k'},
 {'M08_00':'rrho','M12_00':'-rlambda*k'},
 {'M04_00':'rlambda*k^2*q','M12_00':'-rlambda*k'},
 {'M06_00':'rlambda*k/q','M12_00':'-rlambda*k'},
 {'M10_00':'rq0*q','M11_00':'-rq0*q*k'},
 {'M07_00':'rq0*rk*q','M04_00':'-rq0*rk*q*k'}]

source_definitions={
 'k':'sqrt(x)','Q':'1-b','theta':'1-h','H':'h/sig','e':'x*h/sig',
 'K':'1+x','Rs':'1+x*r','r':'b*z','E':'1-e',
 'M':'Rs-x*z*theta','J':'1+x*h-sig*(1-r)',
 'P':'sig+K+sig*z*(1+b*x+2*x*h)+sig*b*x^2*h*z^2+x^2*h*theta*z',
 'p':'(z*(1+x*h)*E+(1-r)*(1+x*h*z))/(2*P)',
 'l':'(1+sig*r)*(e*Rs*(1+sig*z)+M*E)/(2*P)',
 'kappa':'K*(1+sig*r)*(1+sig*z)/(2*M*P)',
 'D':'(K*(1+sig*z)*(1+sig*r)-J*P)/(2*J*P)',
 'beta':'(K*(1+sig*r)*(1+x*h*z)-theta*P)/(2*theta*P)',
 'm':'kappa*x*z*e','C':'K*z*(1+sig*r)*E/(2*(1-r)*P)',
 'a':'p/b','ell':'l/Q','Aevent':'p/(b*z)','alpha':'l/(z*Q)',
 'nu':'(2*D-2*beta*theta+h)/sig',
 'lam':'(8*D*h+4*beta*theta^2+2*h^2)/sig-4*b*p*sig*z',
 'KA':'4*kappa*H-lam','dG':'4*D/sig^2','rho':'lam*x/sig',
 'S':'dG+rho','L0':'lam*sig','q0':'lam*k','c':'max(0,1-KA)',
 'a0':'k/sig','rdG':'sqrt(dG)','rrho':'sqrt(rho)','rsig':'sqrt(sig)',
 'rlambda':'sqrt(lam)','rq0':'sqrt(q0)','rk':'sqrt(k)','q':'sqrt(z)',
 'gs':'h-Q','R0':'-4*(2*(C+D)*(r-1)+2*p+r-1)',
 'Rx':'-4*(2*D*b+(1+2*m)*gs-2*beta*theta)','Bx':'8*kappa*e',
 'Ay':'8*(p-l)+4*z*((1+2*m)*gs+2*D*h+2*p*sig)',
 'chi':'8*(p*b+l*Q)-4*h+4*e*h*(2*D+1)-8*D*h+4*z*(4*D*b*h+4*p*b*sig+(1+2*m)*gs^2+2*beta*theta^2)',
 'k1':'2*(1+2*m)*gs^2+8*D*h*b+8*p*b*sig+4*beta*theta^2',
 'eyu':'2*z*(1+2*m)','eyw':'2*z*(1+2*D)',
 'er0':'z*(8*b*D+8*b*p*sig+4*b*sig*(2*C+2*D+1)*(r-1))',
 'eru':'z*k1','erw':'2*z*b^2*(2*D+1)',
 'ec0':'4*z*(sig*(2*C+1)*(r-1)+2*D*(1+sig*(r-1))+2*p*sig)',
 'ecu':'4*z*(2*D*h+2*p*sig+(1+2*m)*gs)','ecw':'4*z*b*(2*D+1)',
 'Au':'1/k','Av':'1/(sqrt(sig)*q)','Az':'sqrt(sig)/(k*q)',
 'M0':'k^2*q','astar':'lam*k','AC':'c*dG/S',
 'qC0':'1/2-AC*(2/lam+k)','pC':'qC0-AC*delta',
 'C2':'4','C1':'8*astar+8+4/astar','C0':'16+8*astar+4/astar',
 's1':'sqrt(lam*sig)','s2':'sqrt(rho)','s3':'sqrt(lam)*x*q',
 's4':'sqrt(lam*x)/q','g0':'sqrt(lam*x)'}

ordinary_roots=[
 ('2*sqrt(beta)',{'M01_00':'1','M05_00':'h'}),
 ('2*sqrt(beta)',{'O01':'1','O02':'h'}),
 ('sqrt(2*(2*D+1))',{'O03':'1','M08_00':'e'}),
 ('2*sqrt(b*p*z)',{'M05_00':'sig','O04':'-1'}),
 ('2*sqrt(a)',{'O05':'1/q','M09_00':'-b*q'}),
 ('2*sqrt(a)',{'O06':'1/q','M08_00':'-b*q'}),
 ('2*sqrt(a)',{'O07':'1/q','O04':'-b*q'}),
 ('2*sqrt(ell)',{'O05':'1/q','O02':'-1/q','M09_00':'Q*q'}),
 ('2*sqrt(C*z)',{'O06':'1/q','M08_00':'-b*q'}),
 ('sqrt(2+4*m)',{'O01':'1','O02':'h'}),
 ('sqrt(4*a+2*z*(2*C+1)-1)',main_coordinates[1]),
 ('sqrt(4*ell-1)',main_coordinates[2]),
 ('sqrt(L0/2)',{'M05_00':'1','M08_00':'-a0'}),
 ('sqrt((KA+c-1)/2)',main_coordinates[4]),
 ('sqrt((KA+c-1)/2)',main_coordinates[5]),
 ('sqrt((KA+c)/2)',{'M09_00':'q','M07_00':'-q*k'}),
 ('sqrt((KA+c)/2)',{'M10_00':'q*sqrt(2*k)','M11_00':'-q*sqrt(2*k)*k'})]
assert len(ordinary_roots)==17
angular_roots=[
 {'A02':'M0*Au*Av*Az','A01':'-M0'},
 {'A03':'sqrt(2)*M0*Av*Az','A04':'-sqrt(2)*M0*Au'},
 {'A05':'M0*(sqrt(astar)+1/sqrt(astar))*Au*Av',
  'A06':'-M0*sqrt(astar)*Au*Az','A07':'-M0*Az/sqrt(astar)'},
 {'A08':'sqrt(Au)*M0*(sqrt(astar)+1/sqrt(astar))*Au*Av',
  'A09':'-sqrt(Au)*M0*sqrt(astar)*Au*Az',
  'A10':'-sqrt(Au)*M0*Az/sqrt(astar)'},
 {'A01':'sqrt(astar)*M0','A03':'-sqrt(astar)*M0*Av*Az'},
 {'A11':'sqrt(Au)*sqrt(astar)*M0',
  'A12':'-sqrt(Au)*sqrt(astar)*M0*Av*Az'},
 {'A13':'sqrt(2*astar*Av*Az)*M0*Au','A14':'-sqrt(2*astar*Av*Az)*M0'},
 {'A15':'sqrt(Au)*sqrt(2*astar*Av*Az)*M0*Au',
  'A16':'-sqrt(Au)*sqrt(2*astar*Av*Az)*M0'},
 {'A17':'M0*sqrt(Au*Av*Az)*sqrt(C2)*Au'},
 {'A18':'M0*sqrt(Au*Av*Az)*sqrt(C1*Au)'},
 {'A19':'M0*sqrt(Au*Av*Az)*sqrt(C0)'}]
assert len(angular_roots)==11
circuit_roots=[{'C01_0':'sqrt(2*g0)*sqrt(s1)',
                'C02_0':'-sqrt(2*g0)*sqrt(s2)'},
 {'C03_0':'sqrt(2*g0)*sqrt(s3)','C04_0':'-sqrt(2*g0)*sqrt(s4)'},
 {'C05_0':'2*sqrt(g0)*(s1*s2)^(1/4)',
  'C06_0':'-2*sqrt(g0)*(s3*s4)^(1/4)'},
 {'C07_0':'2*sqrt(lam*x)'}]

# Every coefficient references actual labelled features, not an implicit span.
allowed={f['id'] for f in features}
for vec in main_coordinates+angular_roots+circuit_roots:
    assert set(vec)<=allowed
for fac,vec in ordinary_roots: assert set(vec)<=allowed
checks['all_root_feature_references_resolve']=True

# Verify the exported source-coordinate maps themselves, not merely the
# independent identities from which they were transcribed.
SX,SY,SZ,SW,Hlast,eh,hh=sp.symbols('SX SY SZ SW Hlast eh hh')
hh12=[SX,SY,SZ,SW,u*SZ,w*SZ,u*SW,u*w*SW,u*u*SW,
      v**3*SW,v*SW,Hlast]
testvals={f'M{i:02d}_00':ff for i,ff in enumerate(hh12,1)}
testvals.update(dict(zip([f'O{i:02d}' for i in range(1,8)],
    [w*SX,u*w*SZ,u*SX,u*u*w*SW,u*SY,w*SY,u*w*SY])))
pars=dict(src);pars.update({'h':hh,'e':eh,'sqrt':sp.sqrt})
def parse_source(expr):
    return sp.sympify(expr.replace('^','**'),locals=pars)
def linear_map(vec,values):
    return sp.expand(sum(parse_source(coef)*values[key] for key,coef in vec.items()))
expected_main=[SX+hh*u*SZ,SY/q-b*q*u*SW,
    SY/q-w*SZ/q+Q*q*u*SW,rdG*(sig*SZ-u*w*SW),
    q*(u*u-k*k)*SW,k*q*(u-k)*SW,
    rlambda*rsig*u*SZ-rlambda*k*Hlast,
    rrho*u*w*SW-rlambda*k*Hlast,
    rlambda*k*k*q*SW-rlambda*k*Hlast,
    rlambda*k*w*SZ/q-rlambda*k*Hlast,
    rq0*q*(v**3-k*v)*SW,rq0*rk*q*(u-k)*SW]
assert all(sp.expand(linear_map(c,testvals)-expected)==0
           for c,expected in zip(main_coordinates,expected_main))
nx,cc,jj,gg,aa,tt=expected_main[:6]
expected_ordinary=[nx,w*nx,u*(SX+eh*w*SW),u*(sig*SZ-u*w*SW),
    u*cc,w*cc,u*w*cc,u*jj,w*cc,w*nx,cc,jj,
    u*(SZ-a0*w*SW),aa,tt,q*u*(u-k)*SW,q*sp.sqrt(2*k)*v*(u-k)*SW]
assert all(sp.expand(linear_map(c,testvals)-expected)==0
           for (_,c),expected in zip(ordinary_roots,expected_ordinary))
Au,Av,Az,Mzero,astarS,Czero,Cone,Ctwo=sp.symbols('Au Av Az Mzero astarS Czero Cone Ctwo')
Za,Rstar,Ind=sp.symbols('Za Rstar Ind')
pars.update({'Au':Au,'Av':Av,'Az':Az,'M0':Mzero,'astar':astarS,
             'C0':Czero,'C1':Cone,'C2':Ctwo})
avals=[SW,u*w*Za,w*Za,u*SW,u*w*SW,u*Za,Za,v**3*w*SW,
       v**3*Za,v*Za,v*SW,v*w*Za,u*Rstar,Rstar,v**3*Rstar,v*Rstar,
       Ind*v**3*Rstar,Ind*u*Rstar,Ind*v*Rstar]
amap=dict(zip([f'A{i:02d}' for i in range(1,20)],avals))
eang=[Mzero*(Au*Av*Az*u*w*Za-SW),
    sp.sqrt(2)*Mzero*(Av*Az*w*Za-Au*u*SW),
    Mzero*((sp.sqrt(astarS)+1/sp.sqrt(astarS))*Au*Av*u*w*SW
      -sp.sqrt(astarS)*Au*Az*u*Za-Az*Za/sp.sqrt(astarS))]
eang.extend([sp.sqrt(Au)*v*eang[2],
             sp.sqrt(astarS)*Mzero*(SW-Av*Az*w*Za)])
eang.append(sp.sqrt(Au)*v*eang[4])
eang.append(sp.sqrt(2*astarS*Av*Az)*Mzero*(Au*u-1)*Rstar)
eang.append(sp.sqrt(Au)*v*eang[6])
eang.extend([Mzero*sp.sqrt(Au*Av*Az)*Ind*sp.sqrt(Ctwo)*Au*v**3*Rstar,
             Mzero*sp.sqrt(Au*Av*Az)*Ind*sp.sqrt(Cone*Au)*u*Rstar,
             Mzero*sp.sqrt(Au*Av*Az)*Ind*sp.sqrt(Czero)*v*Rstar])
assert all(sp.expand(linear_map(c,amap)-expected)==0
           for c,expected in zip(angular_roots,eang))
gzero,ss1,ss2,ss3,ss4,lams,xs=sp.symbols('gzero ss1 ss2 ss3 ss4 lams xs')
pars.update({'g0':gzero,'s1':ss1,'s2':ss2,'s3':ss3,'s4':ss4,
             'lam':lams,'x':xs})
cf=sp.symbols('cf1:8')
cmap=dict(zip([f'C{i:02d}_0' for i in range(1,8)],cf))
ecircuit=[sp.sqrt(2*gzero)*(sp.sqrt(ss1)*cf[0]-sp.sqrt(ss2)*cf[1]),
    sp.sqrt(2*gzero)*(sp.sqrt(ss3)*cf[2]-sp.sqrt(ss4)*cf[3]),
    2*sp.sqrt(gzero)*((ss1*ss2)**sp.Rational(1,4)*cf[4]
                     -(ss3*ss4)**sp.Rational(1,4)*cf[5]),
    2*sp.sqrt(lams*xs)*cf[6]]
assert all(sp.expand(linear_map(c,cmap)-expected)==0
           for c,expected in zip(circuit_roots,ecircuit))
# Check the seven-channel fixed reserve identity with the evaluated gamma.
qconst=sp.Symbol('qconst')
ajs=sp.symbols('a1:8'); bjs=sp.symbols('b1:8')
AAstar=sp.Integer(2)**30; gg7=1/(7*sp.Integer(2)**38)
split=(qconst-sp.Rational(1,256)+sum(
    ((AAstar+aj)*(gg7+bj)+(AAstar-aj)*(gg7-bj))/2
    for aj,bj in zip(ajs,bjs)))
assert sp.expand(split-(qconst+sum(aj*bj for aj,bj in zip(ajs,bjs))))==0
checks['all_exported_source_coordinate_maps_exact']=True
checks['evaluated_seven_channel_reserve_identity_exact']=True

data={
 'scope':{'source':'0<x,sig<=2^-160;0<z<=1;0<b<1;0<h<1;native gaps;D,beta>=0;lam>0',
          'testing':'abs(us),abs(ws)<=delta; all real state amplitudes',
          'H_upper_bound':None,'global_full_set_pencil':False},
 'constants':{'epsilon':'2^-160','delta':'2^-100','c0':'1/256',
              'Astar':'2^30','gamma7':'1/(7*2^38)'},
 'primary_source_variables':['x','sig','z','b','h'],
 'complete_parent':('4*m*(W-w*X+u*Y-u*w*Z)^2'
  '+4*C*(1+w^2)*(u*W-Y)^2+4*Aevent*(1+u^2)*(1+w^2)*Y^2'
  '+(4*D+2)*(u*w*W+u*X-w*Y-Z)^2'
  '+(4*alpha+2)*(1+u^2)*(Y-w*Z)^2'
  '+(4*beta+2)*(1+w^2)*(X+u*Z)^2'
  '-2*(1+u^2)*(1+w^2)*Z^2'),
 'testing_definitions':{'v':'sqrt(u)','Zplus':'max(Z,0)',
  'Zminus':'max(-Z,0)','rstar':'sqrt(w*W*abs(Z))','Inegative':'1 if Z<0 else 0',
  'gauges':'First angle abs/state sign map in the TeX, then common state sign so W>=0; sign(0)=1.'},
 'source_definitions':source_definitions,
 'features':features,'main_root_order':[str(a) for a in roots],
 'main_root_source_coordinates':main_coordinates,
 'main_gram_upper_triangle':{f'{i+1},{j+1}':str(gram[i,j])
                            for i in range(12) for j in range(i,12)},
 'main_monomials':bt,'main_144_directions':directions,'main_scalar_channels':channels,
 'main_square_recipe':{
  'directions':'Each row vector acts on main_root_order; expand through main_root_source_coordinates.',
  'constant_root':'sqrt(channel.source_constant-c0) times direction dot g',
  'plus_root_j':'sqrt((Astar+a_j)/2) times (direction dot g)*sqrt(gamma7+b_j)',
  'minus_root_j':'sqrt((Astar-a_j)/2) times (direction dot g)*sqrt(gamma7-b_j)',
  'j':'1,...,7'},
 'ordinary_17_roots':[{'source_factor':f,'source_coordinates':c} for f,c in ordinary_roots],
 'angular_11_roots':{'common_source_factor':'sqrt(c*dG/S)',
                    'source_coordinates':angular_roots},
 'circuit_8_roots':{'four_base_source_coordinates':circuit_roots,
  'first_copy_source_factor':'sqrt(pC)',
  'second_copy_source_factor':'sqrt(AC)',
  'second_copy_feature_rule':'replace Cxx_0 by Cxx_1'},
 'checks':checks}
dest=HERE/'c1_positive_lambda_local220.json'
dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
(HERE/'c1_retained_union_check_results.json').write_text(
    json.dumps(checks,indent=2)+'\n',encoding='utf-8')
print(json.dumps(checks,indent=2))
print('Wrote',dest.name)
