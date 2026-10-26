#!/usr/bin/env python3
"""Original exact C1 native coquotients and 4-state parent, SymPy only."""
import sympy as s
k,f,g,a,t=s.symbols('k f g a t',positive=True)
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
ka=k*k;eta=f*f;zet=g*g;vc=ka*eta*zet
Lc=1-a*eta;Hc=1+(1-a)*eta;Qc=1-ka*zet+a*vc
R=1+ka*(1-a*eta);T=1+ka-ka*zet*t; M=R-t;J=T-a;e=k*f/g;Delta=g/(k*f)
P4=1+ka+eta*ka*(a-ka-1)-zet*ka*t+vc*(-zet*ka*t**2+2*ka*t+a*(1-t)+t)
Pa=1+ka-a*eta*ka+zet*((1+ka)-(2+ka)*t)+vc*(-zet*t**2+a*(t-1)+t)
Pd=(1+ka)**2-a*eta*ka**2+zet*ka*((ka+1)-(ka+2)*t)+vc*(-zet*ka*t**2+ka*a*(t-1)+ka*t+ka+1)
c2=2*(1+ka)*(1-t*zet)*Lc*(1+vc)*Hc*Qc
cb=2*(1+ka)*Lc*(1+(1-t)*zet)*Hc*Qc
c3=2*(1+ka)*Lc*Hc*Qc*(1-ka*eta+t*vc)
c5=2*(1+ka)*Lc*Hc*Qc*(1+vc)
c4=2*Qc*P4;cal=2*Lc*Hc*Pa;cd=4*Lc*Qc*Pd
COQ=[c2,cb,c3,c5,c4,cal,cd]
dB=2*a*eta*cd/c3
m=2*a*e**2*c2/(M*c3);C=s.Integer(2);A=2*a*eta**2*c4/c3
al=2*a*Delta**2*k**6*eta**2*cal/c3
D=2*a*c5/(J*c3)-dB/2
be=2*a*e**2*cb/(t*c3)-dB/2
NATIVE=(dB,s.Integer(0),s.Integer(0),m,C,A,D,al,be)
CONTACT=(k,k*f*g,g/f,-k*(1-t*g**2),k*g*Lc/f,s.Integer(1))

def phi_state():
    Du=1+u*u;Dv=1+v*v
    bb=X+u*Z; qa=W+u*Y;ra=u*W-Y;r2=v*qa-bb;r5=ra-v*(Z-u*X)
    H2=4*T*c2; Hb=4*T*cb; H3=4*R*c3;H5=4*R*c5
    return (H2*(R/M)*r2**2+Hb*(R/t)*Dv*bb**2+
       Delta**2*H3*(T/a)*Dv*ra**2+Delta**2*H5*(T/J)*r5**2+
       4*R*T*Delta**2*eta**2*c4*Du*Dv*Y**2+
       4*R*T*Delta**4*k**6*eta**2*cal*Du*(v*Y-Z)**2+
       2*R*T*Delta**2*eta*cd*Du*((v*Y-Z)**2-Dv*Z**2))

def stationary_components():
   return u,v,W,X,Y,Z
