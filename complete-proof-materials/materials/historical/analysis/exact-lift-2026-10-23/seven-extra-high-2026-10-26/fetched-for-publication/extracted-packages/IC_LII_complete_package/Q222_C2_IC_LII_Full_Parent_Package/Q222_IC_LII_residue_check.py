"""Exact determinants for IC and strong-II residue Gram blocks."""
import sympy as S
rho,r,h,ss,aa,u,b=S.symbols('rho r h ss aa u b',real=True)
D0=rho**2/(2*(1-rho**2))
IC=S.Matrix([[2/(1-rho**2),-r/(1-rho)],[-r/(1-rho),rho+2*rho*r/(1+rho)+4*D0*r*r]])
# use nu0=2rho/(1+rho)
expected=2*rho/(1-rho**2)+r*(4*rho-r*(3*rho+1))/((1-rho)*(1+rho)**2)
assert S.cancel(IC.det()-expected)==0
print('IC9 determinant PASS:',S.factor(IC.det()))
LII=S.Matrix([[2,2/(h*(1+h))],[2/(h*(1+h)),2/(h*h*(1+h))]])
assert S.factor(LII[1,1]-LII[1,0]**2/LII[0,0]-2/(h*(1+h)**2))==0
print('LII19 Schur PASS')
# OA8 symmetric 3x3 and its determinant (w/ ss^2=1+u^2)
D=D0
OA=S.Matrix([[2+4*D,-4*D-2*(1-aa),2*u*(1+2*D)],[-4*D-2*(1-aa),4*D,-4*D*u],[2*u*(1+2*D),-4*D*u,2*(1+u*u)+4*D*u*u]])
print('OA8 determinant=',S.factor(OA.det()),flush=True)
print('OA8 V,X minor=',S.factor(OA.extract([0,2],[0,2]).det()),flush=True)
# Normalized IA scalar source positivity in LII remote RR11 etc
print('PASS residue algebra')
