import sympy as sp
s,b,h,t,g=sp.symbols('s b h t g')
x=s*b; B=s*s*b; l=1-h
R=1+h*x; w=R-g; e=s*w; r=1-e
K=1+x; T=1+x*r; jR=T-l/t
Sc=1+B*(l+r/s)
Delta= s*(1/s+K+K*B+2*x*K*r+(B**2*r**2)/s**3+2*B*r/s**2+t*B*(1+1/s+x+B*r*(2-r)/s**2))
Mst=K*(s+r)*(h+t)-s*w*(1+t*B)*l
Qst=1+B*(1+t+h/s+2*r/s+t/s)+B**2*((h+t)/s+(h*(r+1)+t*(1+r-r**2))/s**2)
Kst=(1+s)*(h+t)*K
NN=w*Qst-b*Kst
mg=(1+t*B)*w*T*Sc/(2*Delta)
C=T*Mst/(2*l*Delta)
Ah=t*T*(s+r)*Sc/(2*Delta)
beD=K*(s+r)*(1-t*x*r)*Sc/(2*t*T*Delta)
alsc=(s+r)*Mst/(2*(1-h*B)*Delta)
n=s*s*NN/(2*r*Delta)
n1=s*K*NN/(2*R*Delta)
n2=r*K*K*NN/(2*R*R*Delta)
mu1=mg/R;mu2=mg*g/(R*R);p=w/R
F0=p*p/2-n2-2*n1*p+mu2
AN=p-2*mu1-2*n1
UF=w*w*Delta-(r*K*K+2*s*K*w)*NN+(1+t*B)*w*T*Sc*g
UN=w*(Delta-(1+t*B)*T*Sc)-s*K*NN

def zero(e,label):
    z=sp.cancel(e)
    assert z==0,(label,z)
    print('PASS',label)

zero(Delta-sp.expand(Delta),'Delta polynomial')
print('Delta degree',sp.Poly(sp.expand(Delta),s,b).total_degree())
zero(sp.together(2*Delta*R*R*F0-UF),'F0 numerator')
zero(sp.together(Delta*R*AN-UN),'AN numerator')
for name,poly in [('UF',UF),('UN',UN)]:
    poly=sp.Poly(sp.expand(poly),s,b)
    print(name, 'terms',len(poly.terms()), 'deg',poly.total_degree(), 'b0',sp.factor(poly.eval(b,0)), 's0',sp.factor(poly.eval(s,0)))
    assert poly.eval(b,0)==0
    if name=='UN':assert poly.eval(s,0)==0
for name,q in [('mg',mg),('C',C),('Ah',Ah),('beD',beD),('alsc',alsc)]:
    print('res',name,sp.factor(sp.cancel(q).subs({s:0,g:0})))
# exact degree-two Taylor coefficient at g=0, s,b simultaneously
for name,v in [('C',C),('Ah',Ah),('beD',beD),('n1',n1)]:
    q=sp.factor(v.subs(g,0)); # scale s,b -> tau*s, tau*b and Taylor
    tt=sp.symbols('tt'); z=sp.cancel(q.subs({s:tt*s,b:tt*b}));series=sp.series(z,tt,0,3).removeO();coeff=[sp.factor(series.expand().coeff(tt,i)) for i in range(3)]
    print('TAYLOR',name, coeff)

# Compare canceled table against ORIGINAL native C2 stationary weights.
a=1/s
origDelta=Delta/s
origM=K*(1+a*r)*(h+t)-a*(1+t*B)*e*l
c_a=1+x*r
native = {
  'n': ((1+a*r)*c_a*Sc-r*origDelta)/(2*r*origDelta),
  'mg': a*a*(1+t*B)*e*c_a*Sc/(2*origDelta),
  'C': c_a*origM/(2*l*origDelta),
  'Ah': t*c_a*(1+a*r)*Sc/(2*origDelta),
  'DjR/T': K*(1+a*r)*(1-t*x*r)*Sc/(2*t*jR*origDelta)*(jR/T),
  'alpha': B*B*(1+a*r)*origM/(2*(1-h*B)*origDelta)
}
expected = {'n':n,'mg':mg,'C':C,'Ah':Ah,'DjR/T':beD,'alpha':alsc*s**3*b**2}
for key,expr in native.items():
    zero(expr-expected[key],'true native '+key)
for label,expr in [('mg',mg),('C',C),('A',Ah/h),('D',beD/(jR/T))]:
    z=sp.cancel(expr.subs(g,0))
    print('UNIT-b jet',label,'constant',sp.factor(z.subs(s,0)),
          's-derivative',sp.factor(sp.limit(sp.diff(z,s),s,0)))
