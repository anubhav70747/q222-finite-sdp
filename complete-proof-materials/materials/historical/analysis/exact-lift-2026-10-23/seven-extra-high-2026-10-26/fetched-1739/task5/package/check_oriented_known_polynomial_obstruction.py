import sys,sympy as S
sys.path.insert(0,'/mnt/data/Q222_C1_fixed_admission_2026_10_26')
from c1_native import NATIVE,CONTACT,k,f,g,a,t
vals={k:1,f:1,g:1,a:S.Rational(7,8),t:S.Rational(3,5)}
we=[S.factor(q.subs(vals)) for q in NATIVE]
uc,vc,wc,xc,yc,zc=[S.factor(q.subs(vals)) for q in CONTACT]
d,la,n,m,C,A,D,al,be=we;B=2*C+2*D+d
Sx=S.factor((4*be+2*d)*vc**2*xc**2)
Sy=S.factor(4*A*uc**2*yc**2)
Sz=S.factor(4*be*uc**2*vc**2)
Ub=S.factor(64*m*A*(2*be+d)*B*uc**2*vc**2*wc**4/((2*be+d+2*m)*(2*A+B)))
print('strict actual coefficients:',we)
print('Sx,Sy,Sz=',Sx,Sy,Sz)
print('sqrt side lengths',S.N(S.sqrt(Sx),12),S.N(S.sqrt(Sy),12),S.N(S.sqrt(Sz),12))
print('gamma upper',S.N(S.sqrt(Ub),12),'square limit',S.N(Ub,12))
print('Fourier lower dominant X',S.N((S.sqrt(Sx)-S.sqrt(Sy)-S.sqrt(Sz))**2,12))
print('lower>upper ?', S.N((S.sqrt(Sx)-S.sqrt(Sy)-S.sqrt(Sz))**2 - S.sqrt(Ub),12))
