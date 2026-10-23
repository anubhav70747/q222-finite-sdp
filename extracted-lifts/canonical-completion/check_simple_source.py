"""Full exact certificate for the small valid C1 polynomial-selector obstruction.
Python3 + SymPy. No solver or candidate search; the source is literal data.
"""
from pathlib import Path
import json
import sympy as s

HERE=Path(__file__).resolve().parent
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
m,C,A,D,al,be,d=s.symbols('m C A D alpha beta d')
Du,Dv=1+u*u,1+v*v
F=(4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y*Y
 +(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2*d)*Du*(v*Y-Z)**2
 +(4*be+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z*Z)
variables=(d,m,C,A,D,al,be)
values=tuple(map(s.Rational,('359','75/16','60','416','41/2','831/2','151/2')))
source=dict(zip(variables,values))
uc,vc=s.Rational(1,2),s.Integer(2)
psi=s.Matrix([2,-s.Rational(1,4),s.Rational(1,4),1])
contact={u:uc,v:vc,**dict(zip((W,X,Y,Z),psi))}
checks={}
def ck(name,expr):
    if isinstance(expr,s.MatrixBase):
        ok=all(s.cancel(a)==0 for a in expr)
    else:
        ok=s.cancel(expr)==0
    if not ok:raise AssertionError(name)
    checks[name]=True

assert all(a>0 for a in values)
checks['strictly positive six C1 event weights']=True
P=s.hessian(F,(W,X,Y,Z))/2
Pc=P.subs(source)
Fc=F.subs(source)
ck('complete four-state kernel',Pc.subs({u:uc,v:vc})*psi)
ck('stationary Alice angle',s.diff(Fc,u).subs(contact))
ck('stationary Bob angle',s.diff(Fc,v).subs(contact))
ck('normalization to d=1',F.subs(source)-359*F.subs(dict(zip(variables,[a/359 for a in values]))))
ck('integer homogeneous state',Pc.subs({u:uc,v:vc})*(4*psi))
# Strict m,C,A transverse frame at Z=0; its row determinant is u.
frame=s.Matrix([[v,-1,u*v],[u,0,-1],[0,0,1]])
ck('m C A transverse frame determinant',frame.det()-u)

dd,mm,cc,aa,ee,alpha,beta=values
T=(dd+2*ee,dd+2*(alpha+mm+aa),dd+2*(cc+ee),dd+2*(alpha+beta+mm),
   dd+2*beta,dd+2*(alpha+cc+aa),dd+2*(beta+mm),dd+2*(cc+aa+ee))
q={(0,1):4*alpha*mm*T[4]*T[7],(0,2):4*mm*ee*T[4]*T[5],
   (1,0):4*alpha*aa*T[6]*T[2],
   (1,2):4*(beta*mm*T[0]*T[5]+cc*ee*T[4]*T[1]),
   (2,1):4*cc*aa*T[0]*T[3],(2,2):4*beta*cc*T[0]*T[1]}
assert all(a>0 for a in q.values())
UU,VV=s.symbols('U V')
weights={e:s.factor(a*(uc*uc)**(e[0]-1)*(vc*vc)**(e[1]-1)) for e,a in q.items()}
h=-sum(weights.values())
L=h+sum(a*UU**(i-1)*VV**(j-1) for (i,j),a in q.items())
ck('Laurent value at contact',L.subs({UU:uc*uc,VV:vc*vc}))
ck('Laurent U balance',s.diff(L,UU).subs({UU:uc*uc,VV:vc*vc}))
ck('Laurent V balance',s.diff(L,VV).subs({UU:uc*uc,VV:vc*vc}))
xx,yy=s.symbols('x y',positive=True)
balanced=1800*(3757645*(xx+1/xx-2)+14022297*(yy+1/yy-2)
              +3732435*(yy/xx+xx+1/yy-3)
              +8179670*(xx*yy+1/xx+1/yy-3))
ck('four elementary balanced circuits',L.subs({UU:xx/4,VV:4*yy})-balanced)
det=Pc.det(method='domain-ge')
ck('entire four-state Laurent determinant',det-16*Du**2*Dv**2*u*u*v*v*L.subs({UU:u*u,VV:v*v}))

# Canonical completion forced by the full contact, with source scaled by359.
pp,rr,ss=s.symbols('p r s')
PP=s.Matrix([[4*mm*v*v,-4*mm*v],[-4*mm*v,4*mm+(4*beta+2*dd)*(1+v*v)]])
NN=s.Matrix([[4*cc*(1+v*v)+4*ee+2*dd,(4*ee+2*dd)*v],[(4*ee+2*dd)*v,(4*ee+2*dd)*v*v]])
JJ=s.Matrix([[pp*v*v,0],[rr*v,ss*v*v]])
av=s.Matrix([17,2]);bv=s.Matrix([12,-18]);Uc=uc*uc
rhs1=-PP.subs(v,vc)*av/Uc;rhs2=-Uc*NN.subs(v,vc)*bv
assert rhs1==s.Matrix([-4800,-38400]) and rhs2==s.Matrix([1200,9600])
eq=list(JJ.subs(v,vc)*bv-rhs1)+list(JJ.subs(v,vc).T*av-rhs2)
forced={pp:-100,rr:2000,ss:1200}
ck('every forced-selector contact equation',s.Matrix(eq).subs(forced))
coef,rhs=s.linear_eq_to_matrix(eq,(pp,rr,ss))
assert coef.rank()==coef.row_join(rhs).rank()==3
checks['rank-three forcing system']=True
G=s.BlockMatrix([[PP,JJ],[JJ.T,NN]]).as_explicit().subs(forced)
minor3=G.subs(v,1).extract([1,2,3],[1,2,3])
assert minor3==s.Matrix([[s.Rational(8235,4),2000,1200],[2000,1280,800],[1200,800,800]])
assert minor3.det()==-412640000
checks['literal negative three-dimensional principal minor']=True
minor2=G.subs(v,1).extract([1,2],[1,2])
assert minor2.det()==-1364800
checks['smaller negative two-dimensional principal minor']=True
negative=s.Matrix([0,1,-1,0])
assert (negative.T*G.subs(v,1)*negative)[0]==-s.Rational(2645,4)
checks['two-coordinate negative quadratic witness']=True
ck('source/selector scaling of the three-minor',minor3.det()/359**3+s.Rational(412640000,359**3))

def entries(M):return [[str(a) for a in row] for row in M.tolist()]
data={'scope':'Explicit valid C1 source with no polynomial canonical Bob selector; no obstruction to algebraic selectors or an extra fixed multiplier.',
      'source_order':['d','m','C','A','D','alpha','beta'],'source':list(map(str,values)),
      'd1_source':[str(a/359) for a in values],
      'contact':{'u':str(uc),'v':str(vc),'state':list(map(str,psi)),'equivalent_integer_state':['8','-1','1','4']},
      'contact_state_matrix':entries(Pc.subs({u:uc,v:vc})),
      'SG6_t':list(map(str,T)),'Laurent_central':str(h),
      'Laurent_outer':{f'{i},{j}':str(a) for (i,j),a in q.items()},
      'balanced_weights_at_contact':{f'{i-1},{j-1}':str(a) for (i,j),a in weights.items()},
      'balanced_circuit_identity':{'scale':'1800','x':'4*U','y':'V/4',
         'x+1/x-2':'3757645','y+1/y-2':'14022297','y/x+x+1/y-3':'3732435','x*y+1/x+1/y-3':'8179670'},
      'transformed_contact':{'a':['17','2'],'b':['12','-18'],'U':'1/4','homogeneous_state_scale_relative_to_Z1':'10'},
      'selector_forcing':{'unknown_order':['p','r','s'],'coefficient_matrix':entries(coef),'rhs':list(map(str,rhs)),
          'rank':3,'solution':['-100','2000','1200'],'P_contact':entries(PP.subs(v,vc)),'N_contact':entries(NN.subs(v,vc))},
      'negative_principal_minor':{'v':'1','zero_based_indices':[1,2,3],'matrix':entries(minor3),'determinant':'-412640000'},
      'negative_two_minor':{'v':'1','zero_based_indices':[1,2],'matrix':entries(minor2),'determinant':'-1364800'},
      'negative_vector':{'vector':['0','1','-1','0'],'value':'-2645/4'},
      'exact_checks':len(checks),'checks':checks}
(HERE/'SIMPLE_VALID_SOURCE_AND_SELECTOR.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'exact_checks':len(checks),'scaled_source':data['source'],'contact_state_matrix':data['contact_state_matrix'],
                  'forced_selector':data['selector_forcing']['solution'],'negative_minor3':'-412640000','negative_minor2':'-1364800'},indent=2))
