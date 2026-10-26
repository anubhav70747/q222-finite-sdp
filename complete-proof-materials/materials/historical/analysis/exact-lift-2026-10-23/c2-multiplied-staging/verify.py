"""Exact finite checks for a multiplied native C2 certificate; no SDP solver."""
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sympy as s

u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
d,n,m,C,A,D,al=s.symbols('d n m C A D alpha')
Du,Dv=1+u*u,1+v*v
state=s.Matrix([W,X,Y,Z])
Rm=v*W-X+u*v*Y-u*Z
RD=u*W+u*v*X-Y-v*Z
F=(4*n*Du*Dv*X**2+4*m*Rm**2+4*C*Dv*(u*W-Y)**2
   +4*A*Du*Dv*Y**2+(4*D+2*d)*RD**2
   +(4*al+2*d)*Du*(v*Y-Z)**2+2*d*Dv*(X+u*Z)**2
   -2*d*Du*Dv*Z**2)
P=s.hessian(F,state)/2
H=P[:3,:3]
b=s.Matrix([-u*v*(4*m+4*D+2*d),u*(4*m+2*d-4*D*v*v),
             v*(4*(D-al)-(4*m+4*al+2*d)*u*u)])
zero=lambda z: s.Poly(s.expand(z),*sorted(z.free_symbols,key=str)).is_zero if z.free_symbols else z==0
checks={}
checks['full_parent_cross_column']=all(zero(P[i,3]-b[i]) for i in range(3))

# Store squared row weights, avoiding symbolic square-root simplification.
rows=[]
def add(k,rs):
    rows.extend((k,s.Matrix(r)) for r in rs)
add(2*d,[(u,u*v,-1),(0,0,v),(0,0,u*v),(0,1,0),(0,v,0)])
add(4*n,[(0,1,0),(0,u,0),(0,v,0),(0,u*v,0)])
add(4*m,[(v,-1,u*v)])
add(4*C,[(u,0,-1),(u*v,0,-v)])
add(4*A,[(0,0,1),(0,0,u),(0,0,v),(0,0,u*v)])
add(4*D,[(u,u*v,-1)])
add(4*al,[(0,0,v),(0,0,u*v)])
Hg=sum((k*r*r.T for k,r in rows),s.zeros(3))
checks['nineteen_row_H_gram']=len(rows)==19 and all(zero(t) for t in H-Hg)
adj=H.adjugate()
cross=sum((ki*kj*ri.cross(rj)*ri.cross(rj).T
           for (ki,ri),(kj,rj) in combinations(rows,2)),s.zeros(3))
checks['all_nine_cross_product_adjugate_entries']=all(zero(t) for t in adj-cross)
Hcore=sum((k*r*r.T for k,r in rows[:5]),s.zeros(3))
checks['core_determinant']=zero(Hcore.det()-8*d**3*u*u*v*v*Du*Dv)

# Verify fraction-free Schur identity with independent generic matrix entries.
h11,h12,h13,h22,h23,h33,b1,b2,b3,g=s.symbols('h11 h12 h13 h22 h23 h33 b1 b2 b3 g')
Hg=s.Matrix([[h11,h12,h13],[h12,h22,h23],[h13,h23,h33]])
bg=s.Matrix([b1,b2,b3]); xi=s.Matrix([W,X,Y])
Pg=Hg.row_join(bg).col_join(bg.T.row_join(s.Matrix([[g]])))
eta=Hg*xi+bg*Z
checks['fraction_free_Schur_identity']=zero(Hg.det()*(state.T*Pg*state)[0]
    -(eta.T*Hg.adjugate()*eta)[0]-Pg.det()*Z*Z)

exps=[(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0)]
triangles=[(exps[0],exps[3],exps[4]),(exps[1],exps[2],exps[5])]
x,y=s.symbols('x y',positive=True)
mon=lambda e:x**e[0]*y**e[1]
r1,r2,r3,tau=s.symbols('rho1 rho2 rho3 tau')
pair=sum(rho*(mon(e)+mon((-e[0],-e[1]))-2)
         for rho,e in zip([r1,r2,r3],exps[:3]))
for label,ss,tri in [('plus',[r1+tau,r2,r3,r3+tau,r2+tau,r1],triangles[0]),
                     ('minus',[r1,r2+tau,r3+tau,r3,r2,r1+tau],triangles[1])]:
    original=sum(si*(mon(e)-1) for si,e in zip(ss,exps))
    checks['circuit_coefficients_'+label]=s.expand(original-pair-tau*(sum(mon(e) for e in tri)-3))==0

# At contact Uc=Vc=1, use u=p^3,v=q^3 to clear every root exponent.
# General contact constants are positive source-only powers and restore by scaling.
p,q=s.symbols('p q')
scalar_support=set()
for k,e in enumerate(exps[:3],1):
    root=p**(3+3*e[0])*q**(3+3*e[1])-p**(3-3*e[0])*q**(3-3*e[1])
    target=p**6*q**6*(p**(6*e[0])*q**(6*e[1])+p**(-6*e[0])*q**(-6*e[1])-2)
    checks['pair_square_'+str(k)]=s.expand(root*root-target)==0
    scalar_support.update(s.Poly(root,p,q).monoms())
for label,tri in zip(['plus','minus'],triangles):
    polys=[]
    for eh in tri:
        for ei,ej in combinations(tri,2):
            # Exponents are 3+eh+2ei, not fractional after this substitution.
            root=(p**(3+eh[0]+2*ei[0])*q**(3+eh[1]+2*ei[1])
                 -p**(3+eh[0]+2*ej[0])*q**(3+eh[1]+2*ej[1]))
            polys.append(root)
            scalar_support.update(s.Poly(root,p,q).monoms())
    target=p**6*q**6*(sum(p**(6*e[0])*q**(6*e[1]) for e in tri)-3)
    checks['nine_square_Hurwitz_'+label]=s.expand(sum(t*t for t in polys)/2-target)==0
checks['scalar_root_third_grid_0_to_6']=all(0<=i<=6 and 0<=j<=6 for i,j in scalar_support)
checks['cross_rows_bidegree_at_most_2']=all(s.Poly(z,u,v).degree(u)<=2 and s.Poly(z,u,v).degree(v)<=2
    for (_,ri),(_,rj) in combinations(rows,2) for z in ri.cross(rj) if z!=0)
checks['eta_coefficients_bidegree_at_most_2']=all(s.Poly(z,u,v).degree(u)<=2 and s.Poly(z,u,v).degree(v)<=2 for z in list(H)+list(b) if z!=0)

def native(av,bv,cv,tv,rv):
    z=cv*cv; Ac=1+z*rv; h=bv-av*av*tv
    g=(1+z)*bv-(1-rv+av*av*z*tv); j=1+z*rv-tv
    Delta=((1+z)*z*bv*bv+2*z*(1+z)*bv*rv+(1+z)*(1+av*av*z)*bv
           +z*z*(1-av*av)*rv*rv+2*z*(1+av*av*z)*rv+1+av*av*z)
    N=(1+z)*(bv+rv)*(bv+av*av*(1-tv))-av*av*(1+av*av*z)*(1-rv)*tv
    vals=[1,((bv+rv)*Ac*(Ac+av*av*z*tv)-rv*Delta)/(2*rv*Delta),
          (1+av*av*z)*(1-rv)*Ac*(Ac+av*av*z*tv)/(2*g*Delta),
          Ac*N/(2*av*av*tv*Delta),av*av*Ac*(bv+rv)*(Ac+av*av*z*tv)/(2*h*Delta),
          (1+z)*(bv+rv)*(bv-av*av*z*rv)*(Ac+av*av*z*tv)/(2*av*av*j*Delta),
          z*z*(bv+rv)*N/(2*(1-z*h)*Delta)]
    return dict(zip([d,n,m,C,A,D,al],vals)),[h,g,j,1-z*h,bv-av*av*z*rv],s.Matrix([bv/av,-cv*rv,cv*(bv/av-av*tv),1])
R=s.Rational
example={d:R(1),n:R(7,66),m:R(5,33),C:R(5,22),A:R(10,33),D:R(16,33),al:R(9,22)}
src,gaps,psi=native(R(1),R(5,4),R(1),R(3,4),R(1,4))
checks['counterexample_is_exact_strict_native_source']=src==example and all(t>0 for t in gaps) and all(t>0 for t in src.values())
contact={**example,u:R(1),v:R(1),W:psi[0],X:psi[1],Y:psi[2],Z:1}
checks['counterexample_six_stationary_equations']=all(s.diff(F,t).subs(contact)==0 for t in state) and all(s.diff(F,t).subs(contact)==0 for t in [u,v])
hn=X+u*Z/4; hm=Rm-u*Z; ha=Y-u*Z/2; hd=RD+v*Z/2
remainder=F-4*n*Du*Dv*hn**2-4*m*hm**2-4*A*Du*Dv*ha**2-4*D*hd**2
negative=s.factor(remainder.subs({**example,u:2,v:1,W:0,X:0,Y:0,Z:1}))
checks['four_normal_subtraction_is_minus_31_over_11']=negative==-R(31,11)

seed={d:R(1),n:R(59,320),m:R(73,230),C:R(23,100),A:R(73,200),D:R(73,230),al:R(23,100)}
seed_formula,seed_gaps,_=native(R(1),R(1),R(1),R(5,8),R(1,5))
checks['C2A_source_reconstructed_from_native_formulas']=seed_formula==seed and all(t>0 for t in seed_gaps)
ratios=[]
for src in [example,seed]:
    hsrc=H.subs(src)
    ratios.append(s.factor(hsrc.subs({u:2,v:1}).det()/hsrc.subs({u:1,v:1}).det()))
checks['both_mixed_multiplier_ratios']=ratios==[R(22125,1408),R(460405,31648)] and ratios[0]!=ratios[1]
result={'checks':checks,'all_pass':all(checks.values()),'row_count':len(rows),
        'maximum_squares':19*18//2+12,'feature_count':13*13*4,
        'scalar_root_exponents_in_u_v':[str(Fraction(min(i for i,j in scalar_support),3)),str(Fraction(max(i for i,j in scalar_support),3))],
        'four_normal_remainder':str(negative),'mixed_multiplier_ratios':[str(t) for t in ratios],
        'scope':'SOS for det(H3)*F at independent positive angles; no SOS multiplier cancellation and no affine full lift.'}
out=Path(__file__).with_name('exact_checks.json')
out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
if not result['all_pass']:
    raise SystemExit(1)
