import sympy as sp
k,eta,zeta,a,t=sp.symbols('k eta zeta a t');Vc=k*eta*zeta;K=1+k
Lc=1-a*eta;Hc=1+(1-a)*eta;R=1+k*(1-a*eta);T=1+k-k*zeta*t;Qc=1-k*zeta+a*Vc
P4=1+k+eta*k*(a-k-1)-zeta*k*t+Vc*(-zeta*k*t*t+2*k*t+a*(1-t)+t)
Pa=1+k-a*eta*k+zeta*((1+k)-(2+k)*t)+Vc*(-zeta*t*t+a*(t-1)+t)
Pd=(1+k)**2-a*eta*k*k+zeta*k*((k+1)-(k+2)*t)+Vc*(-zeta*k*t*t+k*a*(t-1)+k*t+k+1)
c2=2*K*(1-t*zeta)*Lc*(1+Vc)*Hc*Qc
cb=2*K*Lc*(1+(1-t)*zeta)*Hc*Qc
c3=2*K*Lc*Hc*Qc*(1-k*eta+t*Vc)
c5=2*K*Lc*Hc*Qc*(1+Vc)
c4=2*Qc*P4
calpha=2*Lc*Hc*Pa
cd=4*Lc*Qc*Pd
balance=[('SD9.1',cb-c2-zeta*c3),('SD9.2',c3-c5+k*eta*cb),('SD9.3',c2-(1-t*zeta)*c5),
('SD10.Y',2*k*c2+2*(1+Vc)*Lc*c4-2*k*k*zeta*Qc*calpha-Qc*cd),
('SD10.S',k*Lc*c2+c5+(1+Vc)*Lc*Lc*c4+k*Qc*Qc*calpha+sp.Rational(1,2)*(k*zeta*Lc*Lc-2*Lc-eta)*cd),
('SD10.V',k*R*c2+t*k*Vc*cb+a*Vc*c3+T*c5+K*Vc*Lc*Lc*c4-k*k*zeta*K*Qc*Lc*calpha-K*sp.Rational(1,2)*(Qc*Lc+eta)*cd)]
for name,ex in balance:
 z=sp.Poly(sp.expand(ex),k,eta,zeta,a,t)
 print(name,'OK' if z.is_zero else 'FAIL',len(z.terms()),flush=True)
