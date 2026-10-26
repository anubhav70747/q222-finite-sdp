import sympy as s
S,k,nu,r,h=s.symbols('S k nu r h');ell=1-nu;theta=1-h;x=k*S;sig=nu*S;L=ell*S;e=k*h/nu;E=1-e
R=1+x*r;J=1+x*h-sig*(1-r);M=R-k*theta/ell;Q=1-L*r;K=1+x
P0=1+S*(nu*ell+2*nu*k*h+k*ell+k*k*h*theta)+S*S*(nu*r*k*ell+nu*r*k*k*h)
P=S*P0;U=K*(1+nu*S*r)
D=(U-J*P0)/(2*J*P0)
beta=(U*(ell+k*h)-theta*P0)/(2*P0*theta)
H0=1+L*(1-r);B=P+x*(J-sig)
Hcal=R*H0+(1-r)*x*h*(2+x*r)-e*(2+x*(1+h+r*theta))+K*e*e
print('sphere beta relation',s.cancel(2*theta*P*beta-h*B+sig*H0)==0,flush=True)
print('sphere total relation',s.cancel(2*P*(J*D+x*theta*beta)-sig**2*Hcal)==0,flush=True)
print('FM D numerator factor S',s.Poly(s.expand(U-J*P0),S).eval(0)==0,flush=True)
