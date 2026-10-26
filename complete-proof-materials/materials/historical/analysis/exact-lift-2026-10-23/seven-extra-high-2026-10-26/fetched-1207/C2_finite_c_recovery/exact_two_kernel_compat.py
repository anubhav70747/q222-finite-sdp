# Recovered checker SOURCE from the previous timed-out run.
# Requires lost first_order_matrices.pkl and two_kernel_grade2.pkl.
# This file was not executed in the recovery session.
import pickle,sympy as s
Z=pickle.load(open('/mnt/data/c2_c_continuation/first_order_matrices.pkl','rb'))
D=pickle.load(open('/mnt/data/c2_c_continuation/two_kernel_grade2.pkl','rb'))
A=D['A'];keys=Z['keys'];
left=A.T.nullspace()
print('left relation count',len(left),flush=True)
e,U,V,W,X,Y,Zz=s.symbols('e U V W X Y Z');a=s.Integer(1);t=s.Rational(3,2);eta=s.Rational(1,4);rho=s.Integer(1)
b=a*a*t+eta*e;zz=1/e;r=rho*e;Ac=1+zz*r;h=b-a*a*t;g=(1+zz)*b-(1-r+a*a*zz*t);j=1+zz*r-t
Dl=(1+zz)*zz*b*b+2*zz*(1+zz)*b*r+(1+zz)*(1+a*a*zz)*b+zz*zz*(1-a*a)*r*r+2*zz*(1+a*a*zz)*r+1+a*a*zz
N=(1+zz)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*zz)*(1-r)*t
n=((b+r)*Ac*(Ac+a*a*zz*t)-r*Dl)/(2*r*Dl)
m=(1+a*a*zz)*(1-r)*Ac*(Ac+a*a*zz*t)/(2*g*Dl)
C=Ac*N/(2*a*a*t*Dl)
A_=a*a*Ac*(b+r)*(Ac+a*a*zz*t)/(2*h*Dl)
D_=(1+zz)*(b+r)*(b-a*a*zz*r)*(Ac+a*a*zz*t)/(2*a*a*j*Dl)
al=zz*zz*(b+r)*N/(2*(1-zz*h)*Dl)
Wc=b/a*W;Xc=-rho*X;Yc=eta/a*Y;Zc=Zz;Vc=a*V;du=U*U+e;dv=Vc*Vc+e
terms=[4*n*du*dv*Xc**2,4*A_*du*dv*Yc**2,4*m*(Vc*Wc+U*Vc*Yc-U*Zc-e*Xc)**2,(4*D_+2)*(U*Wc+U*Vc*Xc-Vc*Zc-e*Yc)**2,(4*al+2)*du*(Vc*Yc-Zc)**2,4*C/e*dv*(U*Wc-e*Yc)**2,2*dv*(2*U*Xc*Zc+e*Xc**2-Zc**2)]
coeff={}
for T in terms:
 p=s.Poly(s.expand(T),U,V,W,X,Y,Zz)
 for powers,co in p.terms():
  k=(powers[0],powers[1],tuple(sorted([('W','X','Y','Z')[i] for i in range(4) for _ in range(powers[i+2])])))
  coeff[k]=coeff.get(k,0)+co
print('all target monomials',len(coeff),'out of map',[(k) for k in coeff if k not in keys],flush=True)
vec=s.Matrix([coeff.get(k,0) for k in keys]+[0]*28)
for ind, L in enumerate(left):
 expr=sum(L[i]*vec[i] for i in range(len(keys)) if L[i]!=0)
 expr=s.factor(s.together(expr))
 print('relation',ind, 'zero' if expr==0 else str(expr)[:450],flush=True)
 if expr!=0:print('FIRST FAILURE at compatibility relation',ind);quit()
print('ALL 13 SYMBOLIC ALL-EPS COMPATIBILITY RELATIONS VANISH',flush=True)
with open('/mnt/data/c2_c_continuation/exact_full_coeffs.pkl','wb') as f:pickle.dump({'coeff':coeff,'vec':vec,'left':left},f)
