"""Fixed testing supports extracted from every full-R1 certificate root.

One finite construction, exact symbolic support and polynomial arithmetic;
no optimizer, degree ladder, random test or received source archive.
All source parameters in root_catalog are coefficients, never functions.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
import json
import sympy as s
HERE=Path(__file__).resolve().parent
u,v,t,r,w,z=s.symbols('u v t r w zeta')
suv,suz,svz,sm=s.symbols('suv suz svz sm')
test=(u,v,t,r,w,z);rad=(suv,suz,svz,sm)
SX,SY,SZ,B,j,T,p,H,alpha,c,cb,sN,sE,tau,rhoQ,sR=s.symbols(
 'SX SY SZ B j T p H alpha c cb sN sE tau rhoQ sR',positive=True)
x,y,fc,qc,rc,gc,kap,Sigma,theta,phi,gamma=s.symbols(
 'x y fc qc rc gc kappa Sigma theta phi gamma',positive=True)
U=SX*u;V=SY*v;Z=SZ*z
Q=U*w-V*Z;Bs=V*w-U*Z;Rs=Z-U*V*w;closed=w-U*V*Z
bc=v*w-u*z/j;Rphys=z-u*v*w/B;Ru=v*(u*u-x)*w
N=t-p*v*w;M=-r-H*v*z
Jalpha=j*alpha*bc;O=s.sqrt(T)*v*(alpha-c*u*u)*Rphys
H1=p*v*v*bc;H2=s.sqrt(T)*p*u*v*bc
PSI=(U*U+V*V)*w*w+(1+U*U*V*V)*Z*Z-4*U*V*w*Z
A27=[Bs,U*Q,U*closed]
Bswap=[Rs,U*Q,U*closed]
E18=[Q,Rs,V*closed]
E=Z*(1-U*V);TT=U*V*(w-Z);dd=U-V;JJ=TT-dd*dd*Z
B29=[E+TT,s.sqrt(2)*dd*w,2*s.sqrt(SX*SY)*suv*(w-Z),
     s.sqrt(2)*JJ,2*s.sqrt(SX*SY)*suv*dd*Z]
psi={1:[(U-V)*w,Z*(1-U*V),s.sqrt(2*SX*SY)*suv*(w-Z)],
    -1:[U*w,V*w,Z,U*V*Z,2*s.sqrt(SX*SY*SZ)*sm]}

G7=s.Matrix([[4,-2,-s.Rational(3,2),0,-s.Rational(1,2),0,s.Rational(1,2)],
 [-2,3,0,-s.Rational(3,2),0,-1,2],[-s.Rational(3,2),0,4,s.Rational(3,2),-2,-3,0],
 [0,-s.Rational(3,2),s.Rational(3,2),4,-s.Rational(3,2),0,-s.Rational(3,2)],
 [-s.Rational(1,2),0,-2,-s.Rational(3,2),3,2,-1],[0,-1,-3,0,2,8,-4],
 [s.Rational(1,2),2,0,-s.Rational(3,2),-1,-4,4]])
G3=s.Matrix([[3,-2,-2],[-2,2,2],[-2,2,3]])
t7=[Z-U*U*V*w,V*w-U*U*V*w,U*w-U*U*V*w,U*Z-U*U*V*w,
 U*V*w-U*U*V*w,U*V*Z-U*U*V*w,U*U*w-U*U*V*w]
t7pi=[V*w-U*U*Z,Z-U*U*Z,U*w-U*U*Z,U*V*w-U*U*Z,
 U*Z-U*U*Z,U*V*Z-U*U*Z,U*U*w-U*U*Z]
Wplus=t7+[2*s.sqrt(SY*SZ)*svz*(1-U)]+[
 s.sqrt(SX*SZ)*suz*f for f in [1-U*V,V-U*V,U-U*V]]
Wpiplus=t7pi+[2*s.sqrt(SY*SZ)*svz*(1-U)]+[
 s.sqrt(SX*SY)*suv*f for f in [w-U*Z,Z-U*Z,U*w-U*Z]]
Wminus=[2*U*w,s.sqrt(3)*V*w,2*U*U*w,s.sqrt(3)*U*U*V*w,s.sqrt(2)*U*V*w,
 2*Z,2*U*Z,s.sqrt(8)*U*V*Z,4*s.sqrt(SX*SY*SZ)*sm,4*U*s.sqrt(SX*SY*SZ)*sm]
Wpiminus=[2*U*w,2*U*U*w,2*V*w,2*U*V*w,s.sqrt(3)*Z,s.sqrt(3)*U*U*Z,
 s.sqrt(2)*U*Z,s.sqrt(8)*U*V*Z,4*s.sqrt(SX*SY*SZ)*sm,4*U*s.sqrt(SX*SY*SZ)*sm]
Wroots={1:Wplus,-1:Wminus};Wpiroots={1:Wpiplus,-1:Wpiminus}
Wpol=4*PSI+4*U*U*(Q*Q+Rs*Rs)-V*V*(U*U-1)**2*w*w
Wp=4*PSI+4*U*U*(Q*Q+Bs*Bs)-Z*Z*(U*U-1)**2
Wmatrix=s.diag(G7,s.ones(1),G3)

def reduce_rad(expr):
    out=0
    for powers,co in s.Poly(s.expand(expr),*rad).terms():
        a,b,c0,d=powers
        assert not any(e%2 for e in powers)
        out+=co*(u*v)**(a//2)*(u*w*z)**(b//2)*(v*w*z)**(c0//2)*(-u*v*w*z)**(d//2)
    return s.expand(out)
checks={}
def check(name,expr):
    value=reduce_rad(expr)
    assert value==0,(name,s.factor(value))
    checks[name]=True
for sign in [1,-1]:
    check('Psi_'+str(sign),sum(f*f for f in psi[sign])-PSI)
    wm=s.Matrix(Wroots[sign]);wpim=s.Matrix(Wpiroots[sign])
    if sign==1:
        check('W19_positive', (wm.T*Wmatrix*wm)[0]-Wpol)
        check('Wpi_positive',(wpim.T*Wmatrix*wpim)[0]-Wp)
    else:
        check('W19_negative',sum(f*f for f in wm)-Wpol)
        check('Wpi_negative',sum(f*f for f in wpim)-Wp)
check('A27',sum(f*f for f in A27)-((1+U*U)*PSI-Rs*Rs))
check('Bswap',sum(f*f for f in Bswap)-((1+U*U)*PSI-Bs*Bs))
check('E18',sum(f*f for f in E18)-((1+V*V)*PSI-V*V*Bs*Bs))
check('B29',sum(f*f for f in B29)-(2*(PSI+V*V*Q*Q+U*U*Bs*Bs)-Rs*Rs))
# Homogeneous RFN20 positive remainder from the swapped W19.
RB=8*U*U*PSI+8*U**4*Q*Q+(8*U**4+2*(U*U-1)**2)*Bs*Bs-V*V*(U*U-1)**2*w*w
check('RFN20',2*U*U*Wp+((U*U-1)*(U*Z-Bs))**2-RB)

# Every exact loss polynomial, before source positive/negative parts.
ss=u*u;vv=v*v
radPsi=(1-phi)*4*x*x*y/fc+phi*(8*theta*x*y*ss/fc+4*(x*x+x*ss+ss*ss)/Sigma)
radQ=(1-phi)*4*x*y*ss/qc+phi*(8*theta*y*ss*ss+4*(1-theta)*x*y*ss)/qc
radE=phi*(8*ss*ss+2*(ss-x)**2)/Sigma
radR=(1-phi)*4*x*y*ss/rc
effectiveRu=16+32*(1-phi)*(alpha/B)**2
Cbar=2*j*T*T*vv*(alpha*alpha+c*c*ss*ss)
hybrid=[phi*((1+theta)/fc+ss/(Sigma*y*x)),2*phi*theta*vv/(qc*y),
        2*phi*ss/(Sigma*y*x),(1-phi)/rc]
loss={
 'Psi':effectiveRu*radPsi+16*phi*j*j*alpha*alpha*y*(1+ss/x)/fc
       +16*Cbar*hybrid[0]+2*kap*hybrid[0]+gc/fc
       +8*gamma*(y*vv+vv*vv+T*ss*(y+vv)),
 'Q':effectiveRu*radQ+16*Cbar*hybrid[1]+2*kap*hybrid[1]+2*kap*vv/(qc*y),
 'E':effectiveRu*radE+16*Cbar*hybrid[2]+2*kap*hybrid[2]+8*gamma*(vv*vv+T*ss*vv),
 'R':effectiveRu*radR+32*(1-phi)*alpha*alpha*j*T*ss/rc+16*Cbar*hybrid[3]+2*kap*hybrid[3]}
loss_support={name:sorted(p for p,c0 in s.Poly(expr,u,v).terms() if p!=(0,0) and c0!=0) for name,expr in loss.items()}
assert all(i%2==k%2==0 for ex in loss_support.values() for i,k in ex)

main=[sN*(N+cb*bc),sE*bc,Ru,Jalpha,H1,M/tau,O,H2,(v*z-T*u*w)/rhoQ,sR*Rphys]
retained=[N,u*N,v*N,v*(N+p*bc),u*M,v*M,u*v*M,V*Q,V*Bs,V*closed,V*Q-Rs,Q,Rs]
origW=w-u*r;origX=t-u*z;origY=r+u*w;origZ=z+u*t
events=[origW+v*origX,origX,v*origX,v*w-t,u*(v*w-t),r,u*r,v*r,u*v*r,
        origY,v*origY,r+v*z,u*(r+v*z),v*origY-origZ,t,u*t,v*t,u*v*t]
radclass={(0,0,0,0):'plain',(1,0,0,0):'uv',(0,1,0,0):'uwz',
          (0,0,1,0):'vwz',(0,0,0,1):'negative_uv_wz'}
catalog=[];base={};gates={};matrix={}
def record(sign,family,expr):
    exps=defaultdict(set)
    for mon,co in s.Poly(s.expand(expr),*(test+rad)).terms():
        if co==0:continue
        cls=radclass.get(mon[6:]);assert cls is not None,(family,mon)
        assert not (sign<0 and cls in ['uwz','vwz'])
        assert not (sign>0 and cls=='negative_uv_wz')
        exps[cls].add(mon[:6])
        base.setdefault((sign,cls),set()).add(mon[:6])
    assert len(exps)==1,(family,exps)
    cls=next(iter(exps))
    catalog.append({'chart':sign,'family':family,'root':str(s.expand(expr)),
                    'radical_class':cls,'monomials':[list(m) for m in sorted(exps[cls])]})
    return cls,exps[cls]
for sign in [1,-1]:
    for name,roots in [('core',main+psi[sign]),('retained',retained),('event',events),
                       ('radial_W',Wroots[sign]),('radial_U_A27',[U*f for f in A27]),
                       ('radial_U_Wpi',[U*f for f in Wpiroots[sign]]),
                       ('radial_RB_square',[(U*U-1)*(U*Z-Bs)]),
                       ('Jalpha_difference',[alpha*Ru/B-alpha*u*Rphys]),('Jalpha_Bswap',Bswap),
                       ('O_v_A27',[v*f for f in A27]),('O_u2v_A27',[u*u*v*f for f in A27]),
                       ('O_v_B29',[v*f for f in B29]),('O_u2v_B29',[u*u*v*f for f in B29]),
                       ('O_majorization',[s.sqrt(T)*v*(alpha+c*u*u)*Rphys]),
                       ('signed_A27',A27),('signed_B29',B29),
                       ('helper_v_E18',[v*f for f in E18]),('helper_u_E18',[u*f for f in E18])]:
        for i,f in enumerate(roots):record(sign,name+'_'+str(i+1),f)
    pivots={'Psi':psi[sign],'Q':[Q],'E':[bc],'R':[Rs]}
    for name,roots in pivots.items():
        for i,f in enumerate(roots):
            cls,support=record(sign,'residual_'+name+'_'+str(i+1),f)
            for a,b in loss_support[name]:
                gates.setdefault((sign,(a,b),cls),set()).update(support)
                # All negative-coefficient branches are included, even when
                # additional source sign information could delete some.
                record(sign,'residual_shift_'+name+'_'+str((a,b,i)),u**(a//2)*v**(b//2)*f)
    for ex,indices in [((1,1),[0,1,5,8]),((0,2),[5,8]),((2,2),[5,8])]:
        support=set()
        for index in indices:
            poly=s.Poly(s.expand(main[index]),*test)
            support.update(m for m,c0 in poly.terms() if c0!=0)
        matrix[sign,ex]=support

def mono(*indices):
    m=[0]*6
    for i in indices:m[i]+=1
    return tuple(m)
zero=(0,)*6
weights={'plain':{zero:F(1)},'uv':{mono(0,1):F(1)},'uwz':{mono(0,4,5):F(1)},
         'vwz':{mono(1,4,5):F(1)},'negative_uv_wz':{mono(0,1,4,5):F(-1)}}
def times_gate(q,ex,sign):
    out=defaultdict(F,q)
    a,b=ex
    for m,c0 in q.items():
        n=list(m);n[0]+=a;n[1]+=b
        out[tuple(n)]+=sign*c0*F(256,9)**a*F(16)**b
    return {m:c0 for m,c0 in out.items() if c0}
blocks=[]
def block(chart,kind,cls,support,q,gate=None):
    blocks.append({'chart':chart,'kind':kind,'radical_class':cls,'gate':gate,
      'monomials':[list(m) for m in sorted(support)],
      'weight':[[list(m),str(c0)] for m,c0 in sorted(q.items())]})
for chart in [1,-1]:
    for (sgn,cls),support in sorted(base.items()):
        if sgn==chart:block(chart,'base',cls,support,weights[cls])
    for (sgn,ex,cls),support in sorted(gates.items()):
        if sgn==chart:block(chart,'scalar-gate',cls,support,times_gate(weights[cls],ex,-1),list(ex))
    for (sgn,ex),support in sorted(matrix.items()):
        if sgn==chart:
            for sign in [1,-1]:block(chart,'matrix-gate','plain',support,times_gate(weights['plain'],ex,sign),[list(ex),sign])

maxdeg=[max(m[k] for b in blocks for m in b['monomials']) for k in range(6)]
assert all(a<=b for a,b in zip(maxdeg,[4,3,1,1,1,1])),maxdeg
chartcounts={}
for sign in [1,-1]:
    chartcounts[sign]={kind:sum(len(b['monomials']) for b in blocks if b['chart']==sign and b['kind']==kind)
                       for kind in ['base','scalar-gate','matrix-gate']}
    chartcounts[sign]['total']=sum(chartcounts[sign].values())
grid={'u':[str(F(i,6)*F(9,256)) for i in range(1,6)],'v':[str(F(i,5)*F(1,16)) for i in range(1,5)],
      't':['1/3','2/3'],'r':['1/3','2/3'],'w':['1/3','2/3'],
      'zeta_positive':['1/3','2/3'],'zeta_negative':['-1/3','-2/3'],
      'points_per_chart':320,'maximum_individual_monomial_degrees':maxdeg,
      'evaluation_Gram_formula':'W_b[e,f]=sum_grid q_b(z)*z^(e+f)/(1+v^2)'}
source_definitions={'SX':'1/sqrt(x)','SY':'1/sqrt(y)','SZ':'1/sqrt(j*T)',
 'c':'B*T*C*H','sN':'sqrt(n+mu+1/2)','cb':'n*j/(n+mu+1/2)',
 'sE':'sqrt(n*(mu+1/2)*j^2/(n+mu+1/2))','tau':'sqrt(T*H*(1-H))',
 'rhoQ':'T*sqrt(x/qc)','sR':'sqrt(rc/(j*T))',
 'phi':'p^2/(1-p+p^2)','theta':'Me/Sigma','Sigma':'Me+fc/y','gamma':'p^2/Sigma',
 'radical_symbols':{'suv':'sqrt(u*v)','suz':'sqrt(u*w*zeta)','svz':'sqrt(v*w*zeta)',
                    'sm':'sqrt(-u*v*w*zeta)'},
 'normalization':'Every fixed feature sqrt(q_b)*monomial is divided by sqrt(1+v^2).'}
family={'testing_cap':{'u':'9/256','v':'1/16'},
 'test_variables':list(map(str,test)),'blocks':blocks,'counts':chartcounts,
 'total_PSD_order_with_nine_source_signs':9+sum(len(b['monomials']) for b in blocks),
 'gram_variables':sum(len(b['monomials'])*(len(b['monomials'])+1)//2 for b in blocks),
 'closure_grid':grid,'source_parameter_definitions':source_definitions}
(HERE/'fixed_family.json').write_text(json.dumps(family,indent=2)+'\n',encoding='utf8',newline='\n')
(HERE/'root_catalog.json').write_text(json.dumps({'source_definitions':source_definitions,'roots':catalog,
 'G7':[[str(q) for q in G7.row(i)] for i in range(7)],'G3':[[str(q) for q in G3.row(i)] for i in range(3)],
 'loss_polynomials':{name:str(s.expand(expr)) for name,expr in loss.items()},
 'loss_support':loss_support},indent=2)+'\n',encoding='utf8',newline='\n')
(HERE/'support_checks.json').write_text(json.dumps({'primitive_polynomial_identities':checks,
 'every_root_has_one_allowed_radical_class':True,'every_root_recorded':len(catalog),
 'maximum_individual_degrees':maxdeg,'counts':chartcounts,
 'matrix_gate_dimensions':{str(k):len(v0) for k,v0 in matrix.items()},
 'loss_support':loss_support,'no_claim_of_299_count':True},indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps({'counts':chartcounts,'total_order':family['total_PSD_order_with_nine_source_signs'],
 'gram_variables':family['gram_variables'],'root_records':len(catalog),'max_degree':maxdeg},indent=2))
