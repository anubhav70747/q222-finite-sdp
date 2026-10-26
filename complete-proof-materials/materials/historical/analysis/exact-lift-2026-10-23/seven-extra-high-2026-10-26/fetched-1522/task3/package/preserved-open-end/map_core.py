"""Deterministic exact 14+2 Gram map, native full parent, no SDP or sample. """
import sympy as S
U,V,W,X,Y,Z=S.symbols('U V W X Y Z'); eps=S.symbols('epsilon', positive=True)
a,t,eta,rho=S.symbols('a t eta rho', positive=True)
n,m,C,A,D,al=S.symbols('n m C A D alpha')
STATES=(W,X,Y,Z)
# Pure source-normalized state: (W_orig,X_orig,Y_orig,Z_orig)=(b/a W,-rho X,eta/a Y,Z), V_orig=aV.
def native_parent(n,m,C,A,D,al):
 b=a*a*t+eta*eps
 vv=a*V; ww=b/a*W; xx=-rho*X; yy=eta/a*Y; zz=Z
 du=U*U+eps; dv=vv*vv+eps
 return (4*n*du*dv*xx**2+4*m*(vv*ww+U*vv*yy-U*zz-eps*xx)**2
  +4*C/eps*dv*(U*ww-eps*yy)**2+4*A*du*dv*yy**2
  +(4*D+2)*(U*ww+U*vv*xx-vv*zz-eps*yy)**2
  +(4*al+2)*du*(vv*yy-zz)**2
  +2*dv*(2*U*xx*zz+eps*xx**2-zz**2))
# source dictionary: set z=c^2=1/eps

def native_weights():
 z=1/eps;b=a*a*t+eta*eps;r=rho*eps
 Ac=1+z*r;h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
 Dl=(1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b+z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z
 Nt=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
 return {n:S.cancel(((b+r)*Ac*(Ac+a*a*z*t)-r*Dl)/(2*r*Dl)),
  m:S.cancel((1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Dl)),
  C:S.cancel(Ac*Nt/(2*a*a*t*Dl)),
  A:S.cancel(a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Dl)),
  D:S.cancel((1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Dl)),
  al:S.cancel(z*z*(b+r)*Nt/(2*(1-z*h)*Dl))}

# feature tensor (angular exponents u,v; state index)
FEATS=[('UW',1,0,0),('UVX',1,1,1),('VZ',0,1,3),('VW',0,1,0),('UZ',1,0,3),('UVY',1,1,2),
 ('UVW',1,1,0),('Z',0,0,3),('UX',1,0,1),('VX',0,1,1),('UY',1,0,2),('VY',0,1,2),('X',0,0,1),('Y',0,0,2)]
ODD=[('sqrtUV W',S.Rational(1,2),S.Rational(1,2),0),('sqrtUV Z',S.Rational(1,2),S.Rational(1,2),3)]

def key(eu,ev,i,j):return (eu,ev,min(i,j),max(i,j))
COLS=[]
for i in range(14):
 for j in range(i,14):COLS.append(('even',i,j))
# Odd kernel: (od0-od1)^2 only to preserve contact; coefficient map includes off-diag terms
COLS.append(('odd',0,1))

def gram_column(col, type_='laplacian'):
 sector,i,j=col
 feat=FEATS if sector=='even' else ODD
 _,u1,v1,s1=feat[i];_,u2,v2,s2=feat[j]
 z=key(u1+u2,v1+v2,s1,s2)
 if type_=='symmetric':return {z:S.Integer(1) if i==j else S.Integer(2)}
 if type_=='laplacian':
  if i==j:return {}
  ki=key(2*u1,2*v1,s1,s1);kj=key(2*u2,2*v2,s2,s2)
  dic={}
  for k,c0 in ((ki,1),(kj,1),(z,-2)):
   dic[k]=dic.get(k,0)+c0
  return {k:S.Integer(c0) for k,c0 in dic.items() if c0}
 raise Exception(type_)

def coefficient_keys():
 # All target monomials up to exponent degree 2, all state pairs; eliminate those impossible by construction.
 return sorted({k for col in COLS for k in gram_column(col,'symmetric')},key=lambda k:(k[0],k[1],k[2],k[3]))
KEYS=coefficient_keys()
KEY_IND={k:i for i,k in enumerate(KEYS)}
def static_matrix(kind='symmetric'):
 return S.MutableSparseMatrix(len(KEYS),len(COLS),{(KEY_IND[k],j):v for j,col in enumerate(COLS) for k,v in gram_column(col,kind).items()})

def target_coeff(nm='symbolic'):
 form=native_parent(n,m,C,A,D,al)
 if nm=='source': form=form.subs(native_weights())
 pp=S.Poly(S.expand(form), U,V,W,X,Y,Z)
 co={}
 for exp,v in pp.terms():
  eu,ev=exp[:2];zz=exp[2:]
  stateindex=[i for i,ii in enumerate(zz) for _ in range(ii)]
  assert len(stateindex)==2
  k=key(eu,ev,*stateindex)
  co[k]=v
 assert set(co).issubset(set(KEYS)), 'off subgrid terms: '+str(set(co)-set(KEYS))
 return S.Matrix([co.get(k,S.Integer(0)) for k in KEYS])
if __name__=='__main__':
 A0=static_matrix('symmetric'); L0=static_matrix('laplacian')
 print('14 even + 2 odd; coefficient rows',len(KEYS),'sym variables',len(COLS),'rank',A0.rank(),'Laplace rank',L0.rank())
 print('nullspace dims',len(COLS)-L0.rank())
