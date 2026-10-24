from sympy import *
a,B,t,h,e = symbols('a B t h e', positive=True)
u,v,W,X,Y,Z = symbols('u v W X Y Z')
r=1-e
# C2 dictionary
ca=1+a*B*r; Sc=1+B*(1-h+a*r); g=1+a*(h*B+r-1); j=ca-(1-h)/t
Dl=a+(1+a*B)+(1+a*B)*B+2*a*B*(1+a*B)*r+a**3*B**2*r**2+2*a**2*B*r+t*B*(1+a+a*B+a**2*B*r*(2-r))
M2=(1+a*B)*(1+a*r)*(h+t)-a*(1+t*B)*e*(1-h)
n2=((1+a*r)*ca*Sc-r*Dl)/(2*r*Dl)
m2=a**2*(1+t*B)*e*ca*Sc/(2*g*Dl)
C2=ca*M2/(2*(1-h)*Dl)
Ah2=t*ca*(1+a*r)*Sc/(2*Dl)
A2=Ah2/h
D2=(1+a*B)*(1+a*r)*(1-a*t*B*r)*Sc/(2*t*j*Dl)
al2=B**2*(1+a*r)*M2/(2*(1-h*B)*Dl)
# C1 native at OP4
kap=a*B; eta=t; zeta=1/a; a1=(1-h)/t; t1=a*e
Vc=kap*eta*zeta; Lc=1-a1*eta; Hc=1+(1-a1)*eta; Qc=1-kap*zeta+a1*Vc
P4=1+kap+eta*kap*(a1-kap-1)-zeta*kap*t1+Vc*(-zeta*kap*t1**2+2*kap*t1+a1*(1-t1)+t1)
Pal=1+kap-a1*eta*kap+zeta*((1+kap)-(2+kap)*t1)+Vc*(-zeta*t1**2+a1*(t1-1)+t1)
Pd=(1+kap)**2-a1*eta*kap**2+zeta*kap*((kap+1)-(kap+2)*t1)+Vc*(-zeta*kap*t1**2+kap*a1*(t1-1)+kap*t1+kap+1)
c2=2*(1+kap)*(1-t1*zeta)*Lc*(1+Vc)*Hc*Qc
cb=2*(1+kap)*Lc*(1+(1-t1)*zeta)*Hc*Qc
c3=2*(1+kap)*Lc*Hc*Qc*(1-kap*eta+t1*Vc)
c5=2*(1+kap)*Lc*Hc*Qc*(1+Vc)
c4=2*Qc*P4; cal=2*Lc*Hc*Pal; cd=4*Lc*Qc*Pd
e1sq=kap*eta/zeta
R1=1+kap*(1-a1*eta); T1=1+kap-kap*zeta*t1; M1=R1-t1; J1=T1-a1
m1=2*a1*e1sq*c2/(M1*c3); A1=2*a1*eta**2*c4/c3
al1=2*a1*(kap**2*zeta/eta)*eta**2*cal/c3
Dtot1=2*a1*c5/(J1*c3); btot1=2*a1*e1sq*cb/(t1*c3); dB=2*a1*eta*cd/c3
Cc1=2
print('M1-g',simplify(M1-g),'J1-j',simplify(J1-j))
# normalize
m1n,A1n,al1n,Dt1n,bt1n,C1n=[simplify(z/dB) for z in (m1,A1,al1,Dtot1,btot1,Cc1)]
beta1=simplify(bt1n-Rational(1,2)); pi_=simplify(1+2*beta1)
kappa=simplify((1+2*n2)/pi_)
x=a*B; s=1/a; w=a*e; b=a**2*B; K=1+x
Theta=n2/B
Pk=ca**2+s*x*(1+t)+s*x**2*(r+h*e+t)
PD=t*x**2*r**2+x*(t-1+e*(1-h-2*t))-1
Pa=1+s-2*e+x*(r**2+s*(r+h*e))
Qc2=1-h*B
