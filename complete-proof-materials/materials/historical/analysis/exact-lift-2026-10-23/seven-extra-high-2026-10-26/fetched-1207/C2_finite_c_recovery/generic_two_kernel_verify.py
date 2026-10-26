# Historical saved checker source, reconstructed from previous run transcript.
# Requires first_order_matrices.pkl and two_kernel_grade2.pkl, missing in current runtime.
# Not executed in recovery session.
import pickle,sympy as s,time
T0=time.time()
print('begin generic exact source-compatibility audit',flush=True)
z=pickle.load(open('/mnt/data/c2_c_continuation/two_kernel_grade2.pkl','rb'))
zz=pickle.load(open('/mnt/data/c2_c_continuation/first_order_matrices.pkl','rb'))
A=z['A'];keys=zz['keys'];left=A.T.nullspace()
e,a,t,eta,rho=s.symbols('e a t eta rho',nonzero=True)
U,V,W,X,Y,Z=s.symbols('U V W X Y Z')
n,m,C,AA,D,al=s.symbols('n m C A D alpha')
b=a*a*t+eta*e;r=rho*e;cc=1/e;Ac=1+cc*r;h=b-a*a*t;g=(1+cc)*b-(1-r+a*a*cc*t);j=1+cc*r-t
Dl=(1+cc)*cc*b*b+2*cc*(1+cc)*b*r+(1+cc)*(1+a*a*cc)*b+cc*cc*(1-a*a)*r*r+2*cc*(1+a*a*cc)*r+1+a*a*cc
Ns=(1+cc)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*cc)*(1-r)*t
source=[((b+r)*Ac*(Ac+a*a*cc*t)-r*Dl)/(2*r*Dl),(1+a*a*cc)*(1-r)*Ac*(Ac+a*a*cc*t)/(2*g*Dl),Ac*Ns/(2*a*a*t*Dl),a*a*Ac*(b+r)*(Ac+a*a*cc*t)/(2*h*Dl),(1+cc)*(b+r)*(b-a*a*cc*r)*(Ac+a*a*cc*t)/(2*a*a*j*Dl),cc*cc*(b+r)*Ns/(2*(1-cc*h)*Dl)]
Wc=b/a*W;Xc=-rho*X;Yc=eta/a*Y;Zc=Z;Vc=a*V;du=U*U+e;dv=Vc*Vc+e
parts=[4*n*du*dv*Xc**2,4*m*(Vc*Wc+U*Vc*Yc-U*Zc-e*Xc)**2,4*C/e*dv*(U*Wc-e*Yc)**2,4*AA*du*dv*Yc**2,(4*D+2)*(U*Wc+U*Vc*Xc-Vc*Zc-e*Yc)**2,(4*al+2)*du*(Vc*Yc-Zc)**2,2*dv*(2*U*Xc*Zc+e*Xc**2-Zc**2)]
for i,p in enumerate(parts):print('part',i,'expansion...',flush=True)
coeff=[{} for _ in parts]
for ix,p in enumerate(parts):
 pp=s.Poly(s.expand(p),U,V,W,X,Y,Z)
 for powers,c0 in pp.terms():
  key=(powers[0],powers[1],tuple(sorted([('W','X','Y','Z')[i] for i in range(4) for _ in range(powers[i+2])])))
  coeff[ix][key]=c0
print('parts expanded',round(time.time()-T0,1),flush=True)
for ii,L in enumerate(left):
 q=0
 for ix in range(7):
  terms=[]
  for k,key in enumerate(keys):
   if L[k] and key in coeff[ix]:terms.append(L[k]*coeff[ix][key])
  if terms:
   piece=s.cancel(s.together(sum(terms)))
   if ix<6:piece=s.cancel(s.together(piece.subs([ (n,source[0]),(m,source[1]),(C,source[2]),(AA,source[3]),(D,source[4]),(al,source[5]) ])))
   q+=piece
 val=s.cancel(s.together(q))
 print('relation',ii,'zero' if val==0 else 'NOTZERO:'+str(s.factor(val))[:900],'time',round(time.time()-T0,1),flush=True)
 if val!=0:quit()
print('ALL GENERIC FACE PARAMETER RELATIONS EXACT',flush=True)
