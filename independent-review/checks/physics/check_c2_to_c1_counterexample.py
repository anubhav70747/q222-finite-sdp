"""A C2 source whose formal C1 image (map OP4) is not positive near the contact.

Native C2 coordinates (c2_short.tex:275-292): c=1/4, a=v_c/u_c=2, b=3/2, t=1/13, r=1/32.
Both parents vanish with zero gradient at the contact; the formal C1 parent has a
negative Hessian direction there, so it cannot be a sum of squares.
Here U, V denote the half-angle tangents u, v. Uses ../thin-layer/base.py.
"""
import os, sys
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "thin-layer", "base.py")
exec(open(BASE).read())
import sympy as sp
U,V,W_,X_,Y_,Z_=sp.symbols('U V W_ X_ Y_ Z_')
# exact native rational source
c=sp.Rational(1,4); an=sp.Integer(2); bn=sp.Rational(3,2); tn=sp.Rational(1,13); rr=sp.Rational(1,32)
z=c**2
h_=bn-an**2*tn; g_=(1+z)*bn-(1-rr+an**2*z*tn); j_=1+z*rr-tn
Ac=1+z*rr
Dl_=((1+z)*z*bn**2+2*z*(1+z)*bn*rr+(1+z)*(1+an**2*z)*bn+z**2*(1-an**2)*rr**2+2*z*(1+an**2*z)*rr+1+an**2*z)
n_=((bn+rr)*Ac*(Ac+an**2*z*tn)-rr*Dl_)/(2*rr*Dl_)
print('gaps h,g,j,1-zh,b-a^2 z r:',h_,g_,j_,1-z*h_,bn-an**2*z*rr,' n=',n_,float(n_))
assert min(h_,g_,j_,1-z*h_,bn-an**2*z*rr,n_)>0
aw=1/bn; Bw=z*bn; tw=an**2/bn; hw=h_/bn; ew=1-rr
sub={a:aw,B:Bw,t:tw,h:hw,e:ew}
m1,A1,al1,Dt1,bt1,C1,kap_,n2v=[sp.nsimplify(sp.simplify(x.subs(sub))) for x in (m1n,A1n,al1n,Dt1n,bt1n,C1n,kappa,n2)]
print('n2 (dictionary) == n native:',sp.simplify(n2v-n_)==0)
D1=Dt1-sp.Rational(1,2); be1=bt1-sp.Rational(1,2)
print('pi=1+2beta1 =',sp.nsimplify(2*bt1),float(2*bt1),' beta1=',float(be1),' coeffs m1,C1,A1,D1,al1:',[float(x) for x in (m1,C1,A1,D1,al1)])
Du=1+U**2; Dv=1+V**2
M=V*W_-X_+U*V*Y_-U*Z_; J=U*W_+U*V*X_-Y_-V*Z_
F1=4*m1*M**2+4*C1*Dv*(U*W_-Y_)**2+4*A1*Du*Dv*Y_**2+4*D1*J**2+4*al1*Du*(V*Y_-Z_)**2+4*be1*Dv*(X_+U*Z_)**2+2*(J**2+Du*(V*Y_-Z_)**2+Dv*(X_+U*Z_)**2-Du*Dv*Z_**2)
# C2 parent for sanity
Cn=Ac*((1+z)*(bn+rr)*h_+an**2*(g_+j_+rr*tn*(1+an**2*z)))/(2*an**2*tn*Dl_)
mn=(1+an**2*z)*(1-rr)*Ac*(Ac+an**2*z*tn)/(2*g_*Dl_); An=an**2*Ac*(bn+rr)*(Ac+an**2*z*tn)/(2*h_*Dl_)
Dn=(1+z)*(bn+rr)*(bn-an**2*z*rr)*(Ac+an**2*z*tn)/(2*an**2*j_*Dl_); aln=z**2*(bn+rr)*((1+z)*(bn+rr)*h_+an**2*(g_+j_+rr*tn*(1+an**2*z)))/(2*(1-z*h_)*Dl_)
F2=4*n_*Du*Dv*X_**2+4*mn*M**2+4*Cn*Dv*(U*W_-Y_)**2+4*An*Du*Dv*Y_**2+4*Dn*J**2+4*aln*Du*(V*Y_-Z_)**2+2*(J**2+Du*(V*Y_-Z_)**2+Dv*(X_+U*Z_)**2-Du*Dv*Z_**2)
vars_=[U,V,W_,X_,Y_,Z_]
s0=[c,an*c,bn/an,-c*rr,c*(bn/an-an*tn),1]
pt=dict(zip(vars_,s0))
for name,F_ in (('F2',F2),('F1',F1)):
    val=sp.simplify(F_.subs(pt)); grad=[sp.simplify(sp.diff(F_,x).subs(pt)) for x in vars_]
    print(name,'value at contact:',val,' gradient:',grad)
H1=sp.Matrix(6,6,lambda i,k: sp.diff(F1,vars_[i],vars_[k]).subs(pt))
H2=sp.Matrix(6,6,lambda i,k: sp.diff(F2,vars_[i],vars_[k]).subs(pt))
import numpy as np
for nm,H in (('F2',H2),('F1',H1)):
    A=np.array(H.evalf(30).tolist(),dtype=float)
    # congruence-scale
    d=np.sqrt(np.abs(np.diag(A))); S=A/np.outer(d,d)
    w,vec=np.linalg.eigh(S)
    print(nm,'normalized Hessian eigenvalues:',w)
    if nm=='F1':
        vv=vec[:,0]/d
        # rationalize
        wr=sp.Matrix([sp.nsimplify(round(x/np.abs(vv).max(),3)) for x in vv])
        q=(wr.T*H*wr)[0]
        print('rational direction',list(wr),' w^T H1 w =',q,'=',float(q))
        wd=sp.Matrix(wr)
        eps=sp.Rational(1,1000)
        ptn={x:(s0[i]+eps*wr[i]) for i,x in enumerate(vars_)}
        print('F1 at contact + eps*w (eps=1/1000):',float(F1.subs(ptn)))
