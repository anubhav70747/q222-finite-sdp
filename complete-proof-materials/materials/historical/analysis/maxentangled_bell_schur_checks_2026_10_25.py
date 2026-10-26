"""Exact real Bell-basis and ordinary-source identities; no optimization."""
import json
from pathlib import Path
import sympy as s

I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Z=s.diag(1,-1)
jxx,jxz,jzz,ux,uz=s.symbols('Jxx Jxz Jzz ux uz',real=True)
op=(-(jxx+jzz)*s.eye(4)
    +ux*(s.kronecker_product(X,I)-s.kronecker_product(I,X))
    +uz*(s.kronecker_product(Z,I)-s.kronecker_product(I,Z))
    +jxx*s.kronecker_product(X,X)
    +jxz*(s.kronecker_product(X,Z)+s.kronecker_product(Z,X))
    +jzz*s.kronecker_product(Z,Z))
F=s.Matrix([1,0,0,1])/s.sqrt(2)
fm=s.Matrix([1,0,0,-1])/s.sqrt(2)
ps=s.Matrix([0,1,1,0])/s.sqrt(2)
pm=s.Matrix([0,1,-1,0])/s.sqrt(2)
U=s.Matrix.hstack(F,-fm,ps,pm)
K=2*s.Matrix([[-jxx,-jxz,ux],[-jxz,-jzz,uz],[ux,uz,-jxx-jzz]])
assert U.T*U==s.eye(4)
assert (U.T*op*U-s.diag(s.Integer(0),K)).applyfunc(s.expand)==s.zeros(4)

qx,qz,mxx,mxz,mzx,mzz=s.symbols('qx qz Mxx Mxz Mzx Mzz',real=True)
derivative=(qx*s.kronecker_product(X,I)+qz*s.kronecker_product(Z,I)
            +mxx*s.kronecker_product(X,X)+mxz*s.kronecker_product(X,Z)
            +mzx*s.kronecker_product(Z,X)+mzz*s.kronecker_product(Z,Z))
# q is the combined local-vector derivative; put all of it on Alice,
# since either-party Pauli actions on Phi+ agree in the X,Z plane.
action=U.T*derivative*F
assert (action-s.Matrix([mxx+mzz,-qz,qx,mzx-mxz])).applyfunc(s.expand)==s.zeros(4,1)

d,n,m,C,A,D,al=s.symbols('d n m C A D alpha')
h1=m-C-D; b1=-m+D-al
t01=d/2+al; t10=-d/2-C; t11=-d/2-m-D
quad=t11*(t01+t10)+t01*h1+t10*b1
source_quad=D*(al+d/2)-m*(C+d/2)
assert s.expand(quad+2*source_quad)==0

data={
 'full_real_planar_Bell_basis_identity_exact':True,
 'normal_basis':['-Phi_minus','Psi_plus','Psi_minus'],
 'normal_block':'2*[[-Jxx,-Jxz,ux],[-Jxz,-Jzz,uz],[ux,uz,-Jxx-Jzz]]',
 'derivative_action_normal_components':['-qz','qx','Mzx-Mxz'],
 'ordinary_coefficient_quadric_factor_exact':True,
 'quadric_factor':'ordinary_quadric = -2*(D*(alpha+d/2)-m*(C+d/2))',
 'Hessian_scope':'Necessary 5x5 joint state/angle matrix at one fixed real planar zero; no global sufficiency.',
 'no_solver_or_numerical_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
