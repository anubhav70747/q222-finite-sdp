"""Native IC full-parent normal readback, exact formal polynomial identities."""
import sympy as s
u,w,W,Z,X,Y,z,y,r,n,m,C,A,D,alpha,B,H=s.symbols('u w W Z X Y z y r n m C A D alpha B H',real=True)
Du=1+u*u;Dw=1+w*w;U=u*u;V=w*w
M=W-w*X+u*Y-u*w*Z
J=u*w*W+u*X-w*Y-Z
F=4*n*Du*Dw*X**2+4*m*M**2+4*C*Dw*(u*W-Y)**2+4*A*Du*Dw*Y**2+4*D*J**2+4*alpha*Du*(Y-w*Z)**2+2*(J**2+Du*(Y-w*Z)**2+Dw*(X+u*Z)**2-Du*Dw*Z**2)
graph={X:-r*u*Z,Y:H*w*Z}
LX=s.cancel(s.diff(F,X).subs(graph)/4);LY=s.cancel(s.diff(F,Y).subs(graph)/4)
PX=(2*D+1)*(B/y)-2*m*B/z*V-2*n*r*(1+U)*(1+V)-r*(2*D+1)*U-2*D-r+V*(1-H-r-2*D*H-2*m*(H+r-1))
fY=2*m-2*C*(1+V)-(2*D+1)*V
PY=B*fY+2*A*H*(1+U)*(1+V)+2*C*H*(1+V)+2*D*H*V+2*D*r*U+2*D-2*alpha*(1-H)*(1+U)+2*m*(H+r-1)*U+(H+r-1)*U+H*(V+1)
b=s.symbols('b',positive=True)
R0=(u*W-B*w*Z)/b;R1=(z*W-B*u*w*Z)/b;R2=u*w*W-B/y*Z;T=(U-z)*Z
LXrhs=(2*D+1)*u*R2-2*m*b/z*w*R1+u*Z*PX
LYrhs=b*fY*R0+w*Z*PY
for name,expr in [('IC10-X',LX-LXrhs),('IC10-Y',LY-LYrhs)]:
    check=s.factor(expr.subs(B,b*b)); print(name,str(check)[:300]);assert check==0
# Monic quotients, use polynomial evaluations (no testing denominators)
sval=1/y
FX=-2*n*r*(1+V)-r*(2*D+1)
GX=-2*m*B/z-2*n*r*(1+z)+1-H-r-2*D*H-2*m*(H+r-1)
FY=2*A*H*(1+V)+2*D*r-2*alpha*(1-H)+(2*m+1)*(H+r-1)
GY=-B*(2*C+2*D+1)+2*A*H*(1+z)+2*C*H+2*D*H+H
for name,poly,f,g in [('PX',PX,FX,GX),('PY',PY,FY,GY)]:
    delta=s.factor(poly-(U-z)*f-(V-sval)*g)
    assert s.factor(delta-poly.subs({U:z,V:sval}))==0,(name,delta)
    print(name+' monic split PASS; contact remainder only')
# LII direct state core some root identities separately
print('IC full normal readback PASS')
