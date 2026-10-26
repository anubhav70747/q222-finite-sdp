"""Universal ME-C2 contact obstruction, with a capped polynomial determinant.

Requires SymPy. No numerical optimizer or parameter sweep. A separate
process caps the final reduced 4x4 polynomial determinant at 30 seconds.
"""
import json
from pathlib import Path
import subprocess
import sys
import sympy as s

if len(sys.argv)>1 and sys.argv[1]=='--det-worker':
    from sympy.polys.matrices import DomainMatrix
    a,r,z=s.symbols('a r z')
    raw=json.loads(sys.stdin.read())
    mat=s.Matrix([[s.sympify(v,locals={'a':a,'r':r,'z':z}) for v in row] for row in raw])
    out=DomainMatrix.from_Matrix(mat).det().as_expr()
    print(json.dumps({'determinant':str(s.factor(out))}))
    raise SystemExit

def reduce(seq):
    out=[]
    for x in seq:
        if out and out[-1]==x: out.pop()
        else: out.append(x)
    return tuple(out)
def mul(x,y): return reduce(x[0]+y[0]),reduce(x[1]+y[1])
def prod(p,q):
    out={}
    for x,c in p.items():
        for y,d in q.items():
            w=mul((tuple(reversed(x[0])),tuple(reversed(x[1]))),y)
            out[w]=out.get(w,0)+c*d
    return out

a,r,z=s.symbols('a r z')
c=s.symbols('c',positive=True)
U=1+z; V=1+a*a*z; K=a+r-1+a*r*z
R=1-a*r*z; EE=1+z*(a+r-a*r)
LL=a*a*(1+r*z)+a*r*(1+(2*r-1)*z)+1-r
JJ=(a-1)*r*r*z-2*a*r*z-a-2*r+1
T0=s.Matrix([
 [(r-1)/U,-a*(a+r)/V,0,-(a-1)*K/(U*V),0,-a*(a+r)/V,(r-1)/U],
 [(1+r*z)/U,(a+r)/V,r,(1+a*z)*K/(U*V),r,(a+r)/V,(r-1)/U],
 [(1-r)/U,-a*(a*r*z-1)/V,-r,-(1+a*z)*K/(U*V),r,a*(a*r*z-1)/V,(1+r*z)/U],
 [-(1+r*z)/U,(a*r*z-1)/V,0,-z*(a-1)*K/(U*V),0,-(a*r*z-1)/V,(1+r*z)/U]])
N=s.Matrix([
 [-(a-1)*(a*a*z+2*a*r*z-1)/((a+1)*V),2*a*(a+r)*U*R/((a+1)*V*K),LL/((a+1)*K)],
 [-(a-1)*(1+(2*r-1)*z)/((a+1)*U),-LL/((a+1)*K),-2*(1-r)*V*(1+r*z)/((a+1)*U*K)],
 [-K*EE/(r*U*V),-(a+r)*R/(r*V),(1-r)*(1+r*z)/(r*U)],
 [1,0,0],
 [(a-1)*K*EE/(r*(a+1)*U*V),(a-1)*(a+r)*R/(r*(a+1)*V),JJ/(r*(a+1)*U)],
 [0,1,0],[0,0,1]])

# Reconstruct the contact-action matrix directly from the Pauli operators.
I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Z=s.diag(1,-1)
def obs(t): return ((1-t*t)*Z+2*t*X)/(1+t*t)
A0=s.kronecker_product(Z,I); A1=s.kronecker_product(obs(c),I)
B0=s.kronecker_product(I,Z); B1=s.kronecker_product(I,obs(a*c)); I4=s.eye(4)
hops=[A1-I4,B1-I4,A0*B0-I4,A1*B1-I4,B0-A0,A0*B1-A0,A1*B0-A0]
psi=s.Matrix([1,-r*c,r*c,1])
action=s.Matrix.hstack(*(h*psi for h in hops))
expected=s.diag(2*c*c,2*c,2*c,2)*T0.subs(z,c*c)
assert (action-expected).applyfunc(s.cancel)==s.zeros(4,7)
pivot=T0[:,[0,1,2,4]]*s.diag(U,V,1,1)
assert s.expand(pivot.det()+2*r*r*(a+1)*K)==0
denN=r*(a+1)*U*V*K
Nh=(denN*N).applyfunc(s.cancel)
Th=(U*V*T0).applyfunc(s.cancel)
assert (Th*Nh).applyfunc(s.expand)==s.zeros(4,3)
assert N[[3,5,6],:]==s.eye(3)
print('Pauli action, pivot and explicit kernel checked.',flush=True)

# Reconstruct all six word coefficient matrices, including their factors.
Id=((),()); aw0=((0,),()); aw1=((1,),()); bw0=((),(0,)); bw1=((),(1,))
roots=[{aw1:1,Id:-1},{bw1:1,Id:-1},{mul(aw0,bw0):1,Id:-1},
       {mul(aw1,bw1):1,Id:-1},{bw0:1,aw0:-1},
       {mul(aw0,bw1):1,aw0:-1},{mul(aw1,bw0):1,aw0:-1}]
words=[((),(0,1)),((0,1),()),((0,1),(0,)),((0,1),(0,1)),
       ((0,1),(1,)),((0,1),(1,0))]
labels=['B0B1','A0A1','A0A1B0','A0A1B0B1','A0A1B1','A0A1B1B0']
Cw=[s.Matrix(7,7,lambda i,j,w=w:s.Rational(prod(roots[i],roots[j]).get(w,0)
                     +prod(roots[j],roots[i]).get(w,0),2)) for w in words]
up=[(i,j) for i in range(7) for j in range(i,7)]
qvars=s.symbols('q0:28'); Q=s.zeros(7)
for (i,j),v in zip(up,qvars): Q[i,j]=Q[j,i]=v
table=[Q[1,4]+Q[2,5]+Q[3,6],-Q[0,4]-Q[0,5]-Q[0,6]+Q[2,6]+Q[3,5],
       Q[0,2]-Q[4,6]-Q[5,6]-Q[6,6],Q[2,3],
       Q[0,5]-Q[3,4]-Q[3,5]-Q[3,6],Q[5,6]]
assert all(s.expand(s.trace(w*Q)-eq)==0 for w,eq in zip(Cw,table))

Mrows=[]
for w in Cw:
    b=N.T*w*N
    Mrows.append([b[0,0],b[1,1],b[2,2],2*b[0,1],2*b[0,2],2*b[1,2]])
M=s.Matrix(Mrows)
assert list(M.row(5))==[0,0,0,0,0,1]
assert list(M.row(3))==[N[2,0],0,0,N[2,1],N[2,2],0]
# S23=0, then Q34=0 eliminates S11. These row operations give
# det(M)=-N[2,0]*det(M4), with no normalization ambiguity.
rows=[0,1,2,4]; cols=[1,2,3,4]
M4=s.Matrix([[s.factor(s.cancel(M[i,j]-M[i,0]*M[3,j]/N[2,0]))
              for j in cols] for i in rows])
P4=M4.copy(); factor=s.Integer(1)
row_factors=[]; col_factors=[]
for i in range(4):
    denominator=s.lcm([s.denom(v) for v in P4.row(i)])
    poly=[s.cancel(v*denominator) for v in P4.row(i)]
    common=s.gcd_list(poly)
    P4[i,:]=s.Matrix([[s.cancel(v/common) for v in poly]])
    factor*=common/denominator
    row_factors.append(str(s.factor(common/denominator)))
for j in range(4):
    common=s.gcd_list(list(P4.col(j)))
    P4[:,j]=s.Matrix([s.cancel(v/common) for v in P4.col(j)])
    factor*=common; col_factors.append(str(s.factor(common)))
print('Six word rows exact; reduced polynomial determinant ready.',flush=True)
worker=subprocess.run([sys.executable,__file__,'--det-worker'],
    input=json.dumps([[str(v) for v in row] for row in P4.tolist()]),
    text=True,capture_output=True,timeout=30,check=True)
det4=s.sympify(json.loads(worker.stdout)['determinant'],locals={'a':a,'r':r,'z':z})
expected_det=32*a*a*(a+r)**2*(1-r)*(1+r*z)*R**2*EE/(r**3*(a+1)**3*U*V**3)
assert s.cancel(-N[2,0]*factor*det4/expected_det)==1
data={
 'Pauli_action_T_and_scaled_T0_exact':True,
 'pivot_columns_one_based':[1,2,3,5],
 'scaled_pivot_determinant':'-2*r^2*(a+1)*K/(U*V)',
 'unscaled_pivot_determinant':'-32*r^2*z^2*(a+1)*K/(U*V)',
 'all_twelve_T0N_entries_zero':True,'free_kernel_rows_one_based':[4,6,7],
 'forbidden_word_order':labels,
 'all_six_displayed_word_equations_have_exact_normalization':True,
 'S_coordinate_order':['S11','S22','S33','S12','S13','S23'],
 'determinant_reduction':'det(M)=-N31*det(M4)',
 'reduced_row_factors':row_factors,'reduced_column_factors':col_factors,
 'reduced_polynomial_matrix':[[str(v) for v in row] for row in P4.tolist()],
 'reduced_polynomial_determinant':str(det4),
 'full_determinant_formula_9_exact':True,
 'full_determinant':'32*a^2*(a+r)^2*(1-r)*(1+r*z)*(1-a*r*z)^2*(1+z*(a+r-a*r))/(r^3*(a+1)^3*(1+z)*(1+a^2*z)^3)',
 'strict_domain':'a,z>0;0<r<1; a*r*z<1; K=a+r-1+a*r*z>0',
 'includes_n_zero':True,'final_determinant_worker_timeout_seconds':30,
 'scope':'Full nine-word nonmembership for ME C2 contacts in the strict finite native chart; not singular limits, not general partially entangled C2, and not the appended C1 seed.',
 'no_numerical_optimizer_or_parameter_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print('Universal determinant (9) checked exactly.',flush=True)
