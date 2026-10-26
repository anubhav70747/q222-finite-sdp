import sympy as s
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
Du=1+u**2; Dv=1+v**2
J=u*W+u*v*X-Y-v*Z;B=X+u*Z;R=u*W-Y;V=v*Z
K=J**2+Du*(v*Y-Z)**2+Dv*B**2-Du*Dv*Z**2
CE2=(R+u*v*B)**2+Dv*B**2+Du*v*v*Y**2+Du*u*u*V**2-2*Du*u*V*(W+v*B)
print('CE2',s.expand(K-CE2)==0)
p,w,H=s.symbols('p w H')
CE6=u**2*H**2-2*u**2*H*Y-2*u*H*W+2*u*w*H*X+u*u*w*w*W**2+2*u*u*w*W*X-2*u*w*w*W*Y+(Du+w*w)*(X*X+Y*Y)-2*u*w*X*Y
Krec=s.expand(K.subs({v:1/w,Z:H/w})*w*w)
print('CE6',s.factor(Krec-CE6)==0)
