"""Exact, finite, solver-free checks for the doubly-infinitesimal C1 BN/UN assembly.

The field labels remain separate in the mathematical proof. This checker
verifies the underlying polynomial identities after formal indeterminate
substitution; it does NOT claim to verify asymptotic positivity.
Requires SymPy. No parameter sampling or numerical optimization.
"""
import sympy as s

checks=[]
def zero(expr, label):
    assert s.cancel(expr)==0, label
    checks.append(label)
    print("CHECK",label,flush=True)

# Formal square-root indeterminates: kappa=ks^2, q=qs^2, S=ss^2,
# beta=bs^2, n=ns^2. Thus no mixed-field radical is used.
ks,qs,ss,bs,ns,l,u,w,W,Z,X,tau=s.symbols(
    'ks qs ss bs ns l u w W Z X tau', nonzero=True)
k=ks**2; q=qs**2; S=ss**2; beta=bs**2; n=ns**2
zeta=l/k; x=k*S; delta=w**2-q
N=k*W-l*u*w*Z
G=u*W/ss-l*ss*w*Z
B=ss*delta*Z
Bu=u*delta*Z
N0=W/zeta-u*w*Z
R=w*(u**2-x)*Z
zero(R-(k*ss*G-u*N)/l,'SN7 radial conductor')
zero(delta**2*N**2+2*k*l*delta*w*G*B+2*k*l**2*w**2*B**2-l**2*(w*Bu)**2-(k*delta*W)**2,'SN8 full helper identity')
zero(delta**2*N0**2+2*k*delta/l*G*(w*B)+2*k*(w*B)**2-(w*Bu)**2-(delta*W/zeta)**2,'BN15 pivot budget')

# BN circuit source rows and both exact helper identities.
A1b=bs*u*Z
A2b=bs*W/(zeta*qs)
A3b=bs*(ks*ss/qs)*w*Z
A4b=bs*u*w*W/(zeta*q*ks*ss)
E1b=A1b-A2b; E2b=A3b-A4b
b=beta/q; z=k*beta/q**2; sqrtz=ks*bs/q
zero(w*B+(w+qs)/sqrtz*E2b+w*(w+qs)/l*G,'BN14 free wB')
zero(b*Bu+(w+qs)*(bs/qs)*E1b+b*(w+qs)*N0,'BN17 complete X-normal cancellation')

# UN circuit source rows and every order-of-allocation identity.
A1n=ns*u*Z
A2n=ns*W/(zeta*qs)
A3n=ns*(ks*ss/qs)*w*Z
A4n=ns*u*w*W/(zeta*q*ks*ss)
E1n=A1n-A2n; E2n=A3n-A4n; E3n=A1n-A3n
sqrtzn=ks*ns/q
zero(w*B+(w+qs)/sqrtzn*E2n+w*(w+qs)/l*G,'UN19 free wB')
zero(Bu+(w+qs)*(qs/ns)*E1n+(w+qs)*N0,'UN20 free Bu')
zero(N0-(qs/ns)*(E2n-E1n)+(w-qs)/ns*E3n-(ks/qs)*w*G/l,'UN21 free N0 at theta-infinitesimal')
h=tau*q;e=tau*k/l
Nx=X+h*u*Z; Troot=w*N0+Bu
Gh=-tau*u*Troot+e*w*ss*G
zero(Nx+tau*Troot-X-e*w*W,'UN16 complete normal shift')
zero(Gh**2+u**2*Nx**2+2*tau*u**2*Troot*Nx -u**2*(X+e*w*W)**2
     +2*e*tau*u*w*ss*Troot*G-e**2*S*w**2*G**2,'UN16 positive hard completion')

# The actual C1 physical-source equalities (ES1--ES3), retained with
# independent bounded source variables, prove the pole cancellations.
S,k,nu,h,r=s.symbols('S k nu h r', positive=True)
l=1-nu; q=nu/l; sig=nu*S; x=k*S; L=l*S; theta=1-h; e=k*h/nu
R=1+x*r; K=1+x; H0=1+L*(1-r); J=1+x*h-sig*(1-r)
M=R-theta*k/l
P=sig*(1+L+2*x*h+r*x*L+r*x*x*h)+K*L+x*x*h*theta
P0=s.cancel(P/S); B0=P+x*(J-sig)
Hcal=(R*H0+(1-r)*x*h*(2+x*r)-e*(2+x*(1+h+r*theta))+K*e**2)
beta=s.cancel((h*B0-sig*H0)/(2*theta*P))
D=s.cancel((sig**2*Hcal/(2*P)-x*theta*beta)/J)
m=s.cancel(K*(1+sig*r)*e*S/(2*M*(l/k)*P))
n=s.cancel(beta*theta**2+2*D*h+h**2/s.Integer(2))
A00=s.cancel(m*(l/k)**2-n/q)
Tnu=r-1-nu-nu*S*(1-r)*(r+l)
U=K*(1+sig*r); Mstar=l*M
T=s.expand(l*J*U+J*Mstar*P0-J*Mstar*k*theta*(J-sig)-2*Mstar*U)
T_at0=s.factor(T.subs(k,0))
zero(T_at0-l*nu*S*Tnu,'UN4 polynomial constant term')
Tk=s.cancel((T-T_at0)/k)
assert s.denom(Tk)==1, 'UN4 monic quotient is not polynomial'
checks.append('UN4 polynomial quotient is finite')
A00_repaired=theta*l*H0/(2*P0)+(e*Tk+h*l*S*Tnu)/(2*J*M*P0)
zero(A00-A00_repaired,'UN4 exact radial-mass cancellation')
t0=(D+h/s.Integer(2)-beta*theta)/q
t0_polefree=D/q+(l*H0-k*(h/q)*(J-sig))/(2*P0)
zero(t0-t0_polefree,'UN15 exact 1/q-pole cancellation')
# Native coquotient ratios, independent check of source entry.
a_native=sig*(1-r); t_native=theta*k/l; eta=1/sig; zeta=l/k
V=x*eta*zeta; Qc=1-x*zeta+a_native*V
P4=(1+x+eta*x*(a_native-x-1)-zeta*x*t_native
    +V*(-zeta*x*t_native**2+2*x*t_native+a_native*(1-t_native)+t_native))
zero(P4-((1+x*h)*(1-e)+(1-r)*(x*h+L)),'native C1 fourth coquotient factorization')
c0=(1+x)*(1+sig*r)*(1-e)/(2*P0)
a0=s.cancel(P4/(2*P0))
assert s.cancel(c0.subs({S:0,k:0,nu:0})-s.Rational(1,2))==0
assert s.cancel(a0.subs({S:0,k:0,nu:0})-s.Rational(1,2))==0
checks.append('native c0 and a0 have unit real residues')

print(f'PASS: {len(checks)} exact symbolic identities and polynomial assertions')
for i,label in enumerate(checks,1): print(f'{i:02d}. {label}')
