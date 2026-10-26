"""Exact R2 circuit, current and opposite-normal algebra; formal variables only."""
import sympy as sp
U,V,Z=sp.symbols('U V Z', positive=True)
Q=U-V*Z; B=V-U*Z; R=Z-U*V
Psi=U**2+V**2+Z**2+U**2*V**2*Z**2-4*U*V*Z
E=Z*(1-U*V); T=U*V*(Z-1);delta=U-V
rhs=sp.Rational(3,2)*delta**2+sp.Rational(1,2)*(E-2*T)**2+(3*U*V+3*V**4)*(1-Z)**2+6*Z*(U*V-V**2)**2
assert sp.expand(sp.Rational(3,2)*Psi+3*V**2*(Q**2+B**2)-R**2-rhs)==0
print('SH1 parameter 3/2 exact PASS',flush=True)
assert sp.expand((1+U**2)*Psi-R**2-B**2-U**2*Q**2-U**2*(1-U*V*Z)**2)==0
assert sp.expand((1+V**2)*Psi-Q**2-R**2-V**2*B**2-V**2*(1-U*V*Z)**2)==0
print('SH2 two currents exact PASS',flush=True)
# RL17 inequality from Young 2ab <= a^2/6+6b^2
Rbudget=(sp.Rational(7,4)*Psi+sp.Rational(21,2)*V**2*Q**2+sp.Rational(7,2)*V**2*B**2)-(Z*(V**2-1))**2
expected=sp.Rational(7,6)*(sp.Rational(3,2)*Psi+3*V**2*(Q**2+B**2)-R**2)+(R/sp.sqrt(6)-sp.sqrt(6)*V*Q)**2
assert sp.expand(Rbudget-expected)==0
print('RL17 complete signed Motzkin/opposite square cost exact PASS',flush=True)
# TU17 current identity and all N terms
Un,Bn,cz,nrm,e=sp.symbols('Un Bn cz N e')
cur=(Un-e*Bn-cz*nrm)**2
comp=(Un-Bn-cz*nrm)**2 + 2*(1-e)*Un*Bn -(1-e**2)*Bn**2 -2*(1-e)*cz*Bn*nrm
assert sp.expand(cur-comp)==0
print('TU17 full current completion in unrestricted N exact PASS',flush=True)
# SH4 factorization by compact polynomial identity
j,p,B,H,T,ell=sp.symbols('j p B H T ell')
K=j+1-p; x=j*B; E0=1+B*(H*T+p-1);y=T*B
assert sp.expand(K+H*j*y-j*E0-(1-p)*(1+x))==0
print('SH4 opposite normal source factorization exact PASS',flush=True)
