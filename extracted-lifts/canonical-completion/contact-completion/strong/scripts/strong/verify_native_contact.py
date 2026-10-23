"""Exact native GE1 checks of all six monic-division contact remainders."""
import json,pickle,sympy as S
from pathlib import Path
O=Path(__file__).resolve().parents[2]/'data'
D=pickle.load((O/'lII_generic.pkl').open('rb'))
z,k,H,B,r,cp,dp,g,q,ss,ts,js,al,nu,mu,aa=S.symbols('z k H B r cp dp g q ss ts js al nu mu aa')
y,t=S.symbols('y t')
Ac=1+z*r; Bb=H+y*t; jj=Ac-t; Sc=Ac+y*t
delta=z+(1+z)*Bb+(1+z)*Bb**2+2*z*(1+z)*Bb*r+z**3*r**2+2*z**2*r+y*((1+z)*Bb+z+z**2*r*(2-r))
Gc=Bb+z*(H+r-1)
Mc=(1+z)*(Bb+z*r)*(H+y)-z*(1+y)*(1-r)*y*t
cc=Ac*Mc/(2*delta)
dd=(1+z)*(Bb+z*r)*(Bb-z*y*r)*Sc/(2*delta)
pair=S.factor(cc*jj+dd*t)
source={B:Bb,cp:cc,dp:dd,g:cc*dd/pair,q:dd*t/pair,ss:Ac/Bb,ts:t,js:jj,
al:(Bb+z*r)*Mc/(2*delta),nu:((Bb+z*r)*Ac*Sc-r*delta)/(2*delta),
mu:(1+y)*(1-r)*Ac*Sc/(2*Gc*delta),aa:Ac*(Bb+z*r)*Sc/(2*H*delta)}
receipt=[]
for i,expr in enumerate(D['contact_residuals']):
    print('native residual',i,flush=True)
    E=expr.subs(k**2,y).subs(source,simultaneous=True)
    num,den=S.fraction(S.together(E))
    zero=S.Poly(num,H,z,y,r,t).is_zero
    assert zero,(i,S.factor(num))
    receipt.append({'residual':i,'cleared_numerator_zero':zero})
(O/'NATIVE_CONTACT_VERIFICATION.json').write_text(json.dumps(receipt,indent=2))
print('All six native contact residuals vanish identically.',flush=True)
