"""Exact universal root identities, no source substitutions or numerics."""
import sympy as S
u,w,W,Z,z,B,y,r=S.symbols('u w W Z z B y r',nonzero=True,positive=True)
rho=B/y; U=u*u; V=w*w
q=S.symbols('q',nonzero=True)
# sqrt(B) treated as real positive radical, symbolic
root=S.sqrt(B)
R0=(u*W-B*w*Z)/root
R1=(z*W-B*u*w*Z)/root
R2=u*w*W-rho*Z
T=(U-z)*Z
IC8=B*(U-z)*(V-1/y)*Z**2 - ((U-z)*Z*R2-z*R0**2+u*R0*R1)
IC8a=rho*T*T+(U-z)*T*R2-u*(z*R0-u*R1)*(u*R0-R1)
for name,form in [('IC8',IC8),('IC8a',IC8a)]:
    assert S.factor(form)==0, (name,S.factor(form))
    print(name,'PASS')

v,k,s,W,xi,zz,NN,z,y=S.symbols('v k s W xi zz NN z y',real=True,nonzero=True)
u=S.symbols('u',real=True,nonzero=True)
U=u*u;V=v*v;Qn=u*xi-zz
Rh=(v*Qn -s*u*(V-y)*W)/k
K=u*k*W
R=v*(U-z)*W
T=u*k*(U-z)*W
LII20=T*T-U*R*R -(U-z)/s*Rh*T +u*(U-z)/s*R*Qn
assert S.factor(LII20.subs(y,k*k))==0
print('LII20 PASS')

# LII15 all five, with source relation k^2=y
srel={y:k*k}
ids={
 'LII15a': z*v*(V-y)*W-((u*V*Qn-u*v*k*Rh)/s - (V-y)*R),
 'LII15b':u*v*(V-y)*W-(V*Qn-v*k*Rh)/s,
 'LII15c':u*k*(V-y)*W-(v*k*Qn-y*Rh)/s,
 'LII15d':U*(U-z)*(V-y)*W*W-(-Rh*T+u*R*Qn)/s,
 'LII15e':U*(V-y)**2*W*W-(V*Qn*Qn-2*v*k*Rh*Qn+y*Rh*Rh)/(s*s)
}
for name,form in ids.items():
    x=S.factor(form.subs(srel));assert x==0,(name,x)
    print(name,'PASS')
print('All toric identities pass.')
# WIE5 full-state replacement, no restriction to the graph.
u,w,W,X,Y,Z,z,B,y,r,H=S.symbols('u w W X Y Z z B y r H',nonzero=True,positive=True)
R2=u*w*W-B/y*Z
T=(u*u-z)*Z
E=S.sqrt(B*H)*(w*w-1/y)*Z
J=u*w*W+u*X-w*Y-Z
e=1+z*r-(B-H)/y
aper=S.sqrt(H/B)
# E_X/sqrt N=X+ruZ, E_Y/sqrt L=Y-HwZ
rhs=J+e*Z+r*T+aper*E-u*(X+r*u*Z)+w*(Y-H*w*Z)
assert S.factor(rhs-R2)==0
print('WIE5 exact full-state identity PASS')
