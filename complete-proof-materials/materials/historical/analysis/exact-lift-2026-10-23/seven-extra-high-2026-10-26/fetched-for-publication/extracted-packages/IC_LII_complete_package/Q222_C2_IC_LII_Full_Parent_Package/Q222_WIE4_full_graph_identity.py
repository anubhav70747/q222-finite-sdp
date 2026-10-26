"""WIE4 removal of the entire large D contact-zero source row; exact toric equality."""
import sympy as S
u,w,W,Z,z,B,H,y,r,n,m,C,A,D,alpha,k=S.symbols('u w W Z z B H y r n m C A D alpha k',positive=True)
U=u*u; V=w*w;rho=B/y
R0=(u*W-B*w*Z)/S.sqrt(B);R1=(z*W-B*u*w*Z)/S.sqrt(B)
R2=u*w*W-rho*Z;T=(U-z)*Z;E=S.sqrt(B*H)*(V-1/y)*Z
kk=k;kap=kk-8*D*H*r;a=S.sqrt(H/B)
# scalar graph IC7
G=(4*B*C-kk*z/B)*R0**2+4*B*m/z**2*R1**2+2*(2*C+2*D+1)*R2**2+r*r*(4*n*V+2*(2*D+2*n+1))*T**2 +H/B*(4*A*U+2*(2*A+2*C+2*D+1))*E**2+kk*u/B*R0*R1+(-4*r*(2*D+1)+kk/B)*R2*T-4*a*(2*C+2*D+1)*R2*E
Q=R2-r*T-a*E
Rhs=(4*B*C-kap*z/B)*R0**2+4*B*m/z**2*R1**2+2*(2*C+1)*R2**2+r*r*(4*n*V+4*n+2)*T**2 +H/B*(4*A*(1+U)+4*C+2)*E**2+kap*u/B*R0*R1+(-4*r+kap/B)*R2*T-4*a*(2*C+1)*R2*E
check=S.simplify(S.expand(G-4*D*Q**2-Rhs))
assert check==0, S.factor(check)
print('WIE4 exact large-D square subtraction PASS')
