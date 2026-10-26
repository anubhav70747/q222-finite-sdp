"""Small additional generic interface checks; no search, no file mutations."""
from map_core import *
from two_kernel import M, DIFFCOLS

c=S.symbols('c', positive=True)
u,v,ww,xx,yy,zz=S.symbols('u v ww xx yy zz')
du=1+u*u; dv=1+v*v
F=(4*n*du*dv*xx**2+4*m*(v*ww-xx+u*v*yy-u*zz)**2
   +4*C*dv*(u*ww-yy)**2+4*A*du*dv*yy**2
   +(4*D+2)*(u*ww+u*v*xx-yy-v*zz)**2
   +(4*al+2)*du*(v*yy-zz)**2+2*dv*(xx+u*zz)**2-2*du*dv*zz**2)
b=a*a*t+eta*eps
scaled=F.subs({u:c*U,v:c*a*V,ww:b/a*W,xx:-rho*X/c,yy:eta*Y/(a*c),zz:Z})/c**2
assert S.cancel(scaled.subs(c,1/S.sqrt(eps))-native_parent(n,m,C,A,D,al))==0
print('PASS generic full-parent normalization, with all six independent weights and four states')

weights=native_weights()
assert S.cancel(weights[C].subs(eps,0))==0
assert S.cancel(S.diff(weights[C],eps).subs(eps,0))==0
print('PASS exact C(0)=C_prime(0)=0 required by checker Taylor truncation')
for name,w in weights.items():
    denominator=S.factor(S.denom(w).subs(eps,0))
    assert denominator != 0
    print('SOURCE',name,'denominator_at_zero=',denominator,'limit=',S.factor(w.subs(eps,0)))

for key0,index,scale in [((0,0,1,1),DIFFCOLS.index((5,5)),1),
                         ((0,0,2,2),DIFFCOLS.index((6,6)),1),
                         ((0,0,1,2),DIFFCOLS.index((5,6)),2)]:
    row=M[KEY_IND[key0],:]
    assert [j for j in range(M.cols) if row[j]] == [index]
    assert row[index]==scale
print('PASS exact constant-angular X^2,Y^2,XY rows uniquely force the two tail diagonals and mixed entry')

cols=list(M.rref()[1]); rows=list(M[:,cols].T.rref()[1])
minor=M[rows,cols]; inverse=minor.inv()
assert M[:,cols]*inverse*M[rows,:]==M
print('PASS fixed rational 66-pivot projection is identity on the full coefficient-map image')
