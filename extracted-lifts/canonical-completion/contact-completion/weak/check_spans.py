"""Exact membership of all recovered roots in the fixed34/16 testing spans."""
from pathlib import Path
import sys,runpy,contextlib,io,json
import sympy as s
HERE=Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):g=runpy.run_path(str(HERE/'received_generator.py'))
for name in ['B','T','j','p','H','u','v','t','r','w','zeta','x','y','qz']:
    globals()[name]=g[name]
su,sv,sz=s.symbols('su sv sz',nonnegative=True)
tests=(u,v,t,r,w,zeta,su,sv,sz)
basis=[]
for root,exponents in [(t,[(0,0),(0,1),(1,0),(1,1)]),
                       (r,[(0,0),(0,1),(1,0),(1,1)]),
                       (w,[(0,0),(0,1),(1,0),(0,2),(1,1),(2,0),(0,3),(1,2),(2,1)]),
                       (zeta,[(0,0),(0,1),(1,0),(0,2),(1,1),(1,2),(2,1),(1,3),(2,2)])]:
    basis.extend(u**i*v**k*root for i,k in exponents)
basis.extend([su*sv*w,su*sv*zeta,su*sz,v*su*sz,u*su*sz,u*v*su*sz,sv*sz,u*sv*sz])
exponent_to_index={s.Poly(b,*tests).monoms()[0]:i+1 for i,b in enumerate(basis)}
assert len(exponent_to_index)==34
principal={1,5,10,11,13,15,16,17,18,19,20,22,23,24,27,28}
results={}
def membership(name,expr,require_principal=False):
    expr=expr.xreplace({s.sqrt(u):su,s.sqrt(v):sv})
    poly=s.Poly(s.expand(expr),*tests)
    indices=[];bad=[]
    for monomial,c in poly.terms():
        if c==0:continue
        index=exponent_to_index.get(monomial)
        if index is None:bad.append((list(monomial),str(c)))
        else:indices.append(index)
    found=not bad and (not require_principal or set(indices)<=principal)
    results[name]={'in_fixed_span':found,'indices':sorted(set(indices)),
                   'required_principal':require_principal,'bad_terms':bad}
    assert found,(name,results[name])
for i,root in enumerate(g['f'],1):membership('main_'+str(i),root,True)

N=g['N0'];MM=g['M0'];bc=g['bcur'];q=g['q'];R=g['R'];L=w-u*v*zeta/(B*qz)
retained={
 'mp_N0':N,'u_N0':u*N,'v_N0':v*N,'v_N0_plus_p_b':v*(N+p*bc),
 'u_M0':u*MM,'v_M0':v*MM,'uv_M0':u*v*MM,
 'VQ':v*q,'VB':v*bc,'V_closed_circuit':v*L,
 'VQ_minus_R':-v*q/(T*B)-R,
 'purchase_v_Q':v*q,'purchase_v_R':v*R,'purchase_v_V_closed':v*v*L,
 'purchase_u_Q':u*q,'purchase_u_R':u*R,'purchase_u_V_closed':u*v*L}
for name,root in retained.items():membership(name,root)

# Recover RFN19 only from its actual current printed scientific body.
G7=s.Matrix([[4,-2,-s.Rational(3,2),0,-s.Rational(1,2),0,s.Rational(1,2)],
 [-2,3,0,-s.Rational(3,2),0,-1,2],
 [-s.Rational(3,2),0,4,s.Rational(3,2),-2,-3,0],
 [0,-s.Rational(3,2),s.Rational(3,2),4,-s.Rational(3,2),0,-s.Rational(3,2)],
 [-s.Rational(1,2),0,-2,-s.Rational(3,2),3,2,-1],
 [0,-1,-3,0,2,8,-4],
 [s.Rational(1,2),2,0,-s.Rational(3,2),-1,-4,4]])
G3=s.Matrix([[3,-2,-2],[-2,2,2],[-2,2,3]])
U,V,Zb,ww=s.symbols('U V Zb ww')
t7=s.Matrix([Zb-U*U*V*ww,V*ww-U*U*V*ww,U*ww-U*U*V*ww,
 U*Zb-U*U*V*ww,U*V*ww-U*U*V*ww,U*V*Zb-U*U*V*ww,U*U*ww-U*U*V*ww])
t3=s.Matrix([1-U*V,V-U*V,U-U*V])
Ps=(U*U+V*V)*ww*ww+(1+U*U*V*V)*Zb*Zb-4*U*V*ww*Zb
Qh=U*ww-V*Zb;Rh=Zb-U*V*ww
W19=4*Ps+4*U*U*(Qh*Qh+Rh*Rh)-V*V*(U*U-1)**2*ww*ww
W19sos=(t7.T*G7*t7)[0]+4*V*ww*Zb*(1-U)**2+U*ww*Zb*(t3.T*G3*t3)[0]
assert s.expand(W19-W19sos)==0
minor7=[G7[:i,:i].det() for i in range(1,8)]
minor3=[G3[:i,:i].det() for i in range(1,4)]
assert minor7==list(map(s.Rational,[4,8]))+[s.Rational(101,4),s.Rational(1049,16),s.Rational(907,16),s.Rational(3279,16),s.Rational(743,16)]
assert minor3==[3,2,2]
images={U:u/s.sqrt(x),V:v/s.sqrt(y),Zb:zeta/s.sqrt(qz),ww:w}
for i,root in enumerate(t7,1):membership('RFN19_t7_'+str(i),root.subs(images,simultaneous=True))
# Common factors (y*qz)^(-1/4) or (x*qz)^(-1/4) are pure source.
membership('RFN19_middle_root',2*sv*sz*(1-u/s.sqrt(x)))
for i,root in enumerate(t3,1):
    membership('RFN19_t3_'+str(i),su*sz*root.subs(images,simultaneous=True))
payload={'roots':results,'all_roots_in_fixed_spans':all(x['in_fixed_span'] for x in results.values()),
 'RFN19_homogeneous_identity_exact':True,'RFN19_G7_leading_minors':list(map(str,minor7)),
 'RFN19_G3_leading_minors':list(map(str,minor3)),
 'note':'Constant/source square-root factors do not change these testing spans. No source Gram PSD bounds checked.'}
(HERE/'span_checks.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf8')
print('Exact RFN19 homogeneous SOS and all',len(results),'root memberships passed.')
print('Eleven main-root spans:',[results['main_'+str(i)]['indices'] for i in range(1,12)])
