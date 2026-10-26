import sympy as s
j,B,T,p,H=s.symbols('j B T p H')
x=j*B;y=T*B;K=j+1-p
P=j+1+B*((j+1)**2-p)+T*B*(H*j+j+1+B*(H*j*(2*j+2-p)+p*(2*j+1-p)))+H*j*T**2*B**2*(1+B*p+B*H*j)
L=K+T*((1-H)*j+1-p-H)+B*T*(p*K+H*j*(1-p))+B*H*T**2*((1-H)*j-p)
V=(1+y)*(1+x)*(1+H*y)
U=(1+x)*(1+H*y)*(K+H*j*y)
q=s.Poly(s.expand(P-p*V-B*L-U),j,B,T,p,H)
print('S0 reduced polynomial coefficients',len(q.terms()));assert q.is_zero
print('PASS R34 S0 by short polynomial identity')
E=1+B*(H*T+p-1)
N=(p+B*(-j*j+3*j*p-j-p*p+p)+B*T*(H*j*j+H*j-j*j+j*p-j+p)+B*B*T*(2*H*j*j*p-H*j*j-H*j*p*p+2*H*j*p-2*j*j*p+3*j*p*p-j*p-p**3+p*p)+B*B*H*T*T*j*(H*j-j+p)*(1+B*p))
S2num=-j*p*V+(1+H*y)*(K+H*j*y)*(1-y*((1-H)*j-p))-P*(1-p)+B*(K+H*j*y)*L
q=s.Poly(s.expand(S2num),j,B,T,p,H)
print('S2 reduced polynomial coefficients',len(q.terms())); assert q.is_zero
print('PASS R34 S2 by pure source polynomial identity')
SNnum=(1+y)*(P*p-N)-j*p*V*(1+x)+j*B*B*(K+H*j*y)*L+j*B*P*(-(1-p)+T*(H-1))
q=s.Poly(s.expand(SNnum),j,B,T,p,H)
print('SN reduced polynomial coefficients',len(q.terms())); assert q.is_zero
print('PASS R34 SN by pure source polynomial identity')
