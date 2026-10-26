"""Fast factor-free polynomial certificates for RL (all exact)."""
import sympy as sp
j,p,S,y,ell=sp.symbols('j p S y ell')
B=S*y;H=S*ell;K=j+1-p;x=j*B;W0=(1-H)*j+1-p;E=1+B*(ell+p-1)
P=j+1+B*((j+1)**2-p)+y*(H*j+j+1+B*(H*j*(2*j+2-p)+p*(2*j+1-p)))+H*j*y*y*(1+B*p+B*H*j)
Lstar=S*K+(W0-H)+S*y*(p*K+H*j*(1-p))+S*y*ell*(W0-1)
N=p+B*(-j*j+3*j*p-j-p*p+p)+y*(H*j*j+H*j-j*j+j*p-j+p)+S*y*y*(2*H*j*j*p-H*j*j-H*j*p*p+2*H*j*p-2*j*j*p+3*j*p*p-j*p-p**3+p*p)+H*j*y*y*(H*j-j+p)*(1+B*p)
muN=j*p*(1+y)*(1+x)*(1+H*y)
AN=S*y*y*(K+H*j*y)*Lstar
AC=(1+x)*(1+H*y)*(K+H*j*y)
DN=(1+H*y)*(K+H*j*y)*(1-y*((1-H)*j-p))
# Actual: 2P mu = muN, 2P n = N/(j-p), 2P alpha = AN/E,
# 2P CH = Lstar, 2P A = AC/[S(1-ell)], 2P D = S DN/(1-H)
checks={
'S0': Lstar-S*DN+H*P-AC,
'S1': j*(1+y)*Lstar-j*S*DN-j*(1-H)*P+muN,
'S2': -muN*E*S*y +(Lstar-AC)*E*y +(AN+E*P)*S*y*(ell-1+p)+AN,
'SN': -E*N-E*muN*(1+x)+E*P*p-x*(1-p)*(AN+E*P)+j*AN+y*(E*P*(p-j)+j*H*(AN+E*P)-E*N),
'RL5 p-zero': N.subs(p,0)+y*j*(j+1)*(1+S*(1-ell))*(1+S*ell*y*j/(j+1)),
}
for key,expr in checks.items():
    if key=='RL5 p-zero':expr=sp.together(expr)*(j+1)
    poly=sp.Poly(expr,S,y,j,p,ell)
    ok=poly.is_zero
    print(key,ok,'terms',len(poly.terms()),flush=True)
    if not ok:
        print('FIRST TERMS', poly.terms()[:5]); raise SystemExit(1)
print('RFN21/22 source cancellations and RL5 exact polynomial PASS')
