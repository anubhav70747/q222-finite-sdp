#!/usr/bin/env python3
"""Exact, solver-free verifier and COO exporter for SF1--SF12.

The analytic reverse inclusion uses the degree-bounded matrix Fejer--Riesz
theorem on R.  This program verifies all particular polynomial identities,
rank claims and evaluated affine arrays used in that argument.  It neither
imports an SDP solver nor claims the universal two-angle C2 converse.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
CHECKS = []


def require(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def zero(name, value):
    require(name, s.cancel(value) == 0)


def matrix_zero(name, value):
    require(name, all(s.cancel(x) == 0 for x in value))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')


u, v = s.symbols('u v')
W, X, Y, Z = state = s.symbols('W X Y Z')
d, n, m, C, A, D, alpha = source = s.symbols('d n m C A D alpha')
du, dv = 1+u*u, 1+v*v
P = s.expand(
    n*du*dv*X*X + m*(v*W-X+u*v*Y-u*Z)**2
    + C*dv*(u*W-Y)**2 + A*du*dv*Y*Y
    + (D+d/2)*(u*W+u*v*X-Y-v*Z)**2
    + (alpha+d/2)*du*(v*Y-Z)**2
    + d/2*dv*(X+u*Z)**2-d/2*du*dv*Z*Z)
M = s.hessian(P, state)/2
K = P.diff(d)
zero('complete signed core restricted to WXY', K.subs(Z, 0)
     - ((u*W+u*v*X-Y)**2+du*v*v*Y*Y+dv*X*X)/2)
zero('WXY core row determinant is nonzero at uv nonzero',
     s.Matrix([[u,u*v,-1],[0,0,v],[0,1,0]]).det()+u*v)
for name, sub in [('Alice', {u:-u,Y:-Y,Z:-Z}),
                  ('Bob', {v:-v,X:-X,Z:-Z})]:
    zero('whole parent '+name+' sign gauge',
         P.subs(sub, simultaneous=True)-P)
reciprocal = P.subs({u:1/v,v:1/u,W:Z,Z:W,m:D,D:m,C:alpha,alpha:C},
                   simultaneous=True)
zero('complete reciprocity identity', u*u*v*v*reciprocal-P)

outer = {
    (4,2): A*C*(2*D+d+2*n)*(2*alpha+d+2*m)/4,
    (4,0): A*n*(2*C+2*D+d)*(2*alpha+d+2*m)/4,
    (2,4): C*D*(d+2*n)*(2*A+2*alpha+d+2*m)/4,
    (2,0): A*alpha*(2*C+2*D+d)*(d+2*m+2*n)/4,
    (0,4): D*m*(d+2*n)*(2*A+2*C+2*alpha+d)/4,
    (0,2): alpha*m*(d+2*n)*(2*A+2*C+2*D+d)/4,
}
det_red = s.cancel(M.det(method='domain-ge')/(du*du*dv*dv))
det_terms = s.Poly(det_red,u,v).as_dict()
require('SF3 full seven-term determinant support',
        set(det_terms) == set(outer)|{(2,2)})
for power, value in outer.items():
    zero('SF3 factored determinant coefficient '+str(power),
         det_terms[power]-value)
zero('SF3 central coefficient exact definition',
     det_terms[(2,2)]-M.subs({u:1,v:1}).det()/16+sum(outer.values()))

fixed = {D:m,alpha:C}
a0=A*C*(2*m+d+2*n)*(2*C+d+2*m)/4
b0=C*m*(d+2*n)*(2*A+2*C+d+2*m)/4
for power in [(4,2),(2,0)]:
    zero('SF4 first opposite pair '+str(power), outer[power].subs(fixed)-a0)
for power in [(2,4),(0,2)]:
    zero('SF4 second opposite pair '+str(power), outer[power].subs(fixed)-b0)
rr, tt = s.symbols('r t', positive=True)
U,V=s.symbols('U V', positive=True)
Laurent = sum(val.subs(fixed)*U**(i//2-1)*V**(j//2-1)
              for (i,j),val in det_terms.items())
zero('SF5 exact global angular projection identity',
     Laurent.subs({U:rr*tt,V:rr/tt}, simultaneous=True)
     -Laurent.subs({U:tt,V:1/tt}, simultaneous=True)
     -(rr+1/rr-2)*(a0*tt+b0/tt))
zero('SF5 scalar radial factor is a square over positive r',
     rr+1/rr-2-(rr-1)**2/rr)

H=(u*u*M.subs(fixed).subs(v,1/u)).applyfunc(s.cancel)
require('SF6 entire univariate matrix is polynomial of degree at most four',
        all(s.Poly(e,u).degree()<=4 for e in H if e != 0))
J4=s.diag(1,-1,-1,1)
matrix_zero('whole H signed real-line congruence', H.subs(u,-u)-J4*H*J4)
Swap=s.Matrix([[0,0,0,1],[0,1,0,0],[0,0,1,0],[1,0,0,0]])
matrix_zero('fixed H commutes with W-Z exchange', H-Swap.T*H*Swap)

S,T=s.symbols('S T')
Qstate=s.Matrix([[1,0,0,1],[0,1,0,0],[0,0,1,0],[1,0,0,-1]])
require('SF6 constant state congruence invertible', Qstate.det()==-2)
Cong=(Qstate.T*H*Qstate).applyfunc(s.expand)
R=s.Matrix([
    [(2*C+d+2*m)*u**4+(2*C-d-4*m)*u*u+2*m,
     (d+2*m)*u**3-2*m*u,
     -(2*C+d+2*m)*u**3+(2*m-2*C)*u],
    [(d+2*m)*u**3-2*m*u,
     n*u**4+(2*n+2*m+d)*u*u+n+d/2,
     -(d/2+2*m)*u*u],
    [-(2*C+d+2*m)*u**3+(2*m-2*C)*u,
     -(d/2+2*m)*u*u,
     A*u**4+(2*A+2*C+2*m+d)*u*u+A+2*C+d/2]])
qt=(1+u*u)*((2*C+d+2*m)*u*u+2*m)
matrix_zero('SF6 complete block decomposition including every state normal',
            Cong-s.diag(R,qt))
J3=s.diag(1,-1,-1)
matrix_zero('SF7 parity congruence', R.subs(u,-u)-J3*R*J3)
zero('SF6 separated nonnegative T square coefficient',
     qt-2*m-(2*C+d+4*m)*u*u-(2*C+d+2*m)*u**4)

fplus=s.Matrix([S,u*X,u*Y,u*u*S])
fminus=s.Matrix([X,Y,u*S,u*u*X,u*u*Y])
for i, f in enumerate(fplus):
    zero('SF8 positive parity feature '+str(i),
         f.subs({u:-u,X:-X,Y:-Y},simultaneous=True)-f)
for i, f in enumerate(fminus):
    zero('SF8 negative parity feature '+str(i),
         f.subs({u:-u,X:-X,Y:-Y},simultaneous=True)+f)
all_features=list(fplus)+list(fminus)
require('SF8 complete nine monomial degree-two state-linear space',
        set(all_features)=={z*u**k for z in [S,X,Y] for k in range(3)})

x=s.symbols('x1:10')
x1,x2,x3,x4,x5,x6,x7,x8,x9=x
Qplus=s.Matrix([[2*m,x1,x2,x3],[x1,x4,x5,x6],[x2,x5,x7,x8],
                [x3,x6,x8,2*C+d+2*m]])
l1,l2=-2*m-x1,2*m-2*C-x2
l3,l4=n+m+d/2-x4/2,-d/2-2*m-x5-x9
l5,l6=A+C+m+d/2-x7/2,2*C-d-4*m-2*x3
l7,l8=d+2*m-x6,-2*C-d-2*m-x8
Qminus=s.Matrix([[n+d/2,0,l1,l3,x9],[0,A+2*C+d/2,l2,l4,l5],
                 [l1,l2,l6,l7,l8],[l3,l4,l7,n,0],[x9,l5,l8,0,A]])
target=(s.Matrix([S,X,Y]).T*R*s.Matrix([S,X,Y]))[0]
gram=(fplus.T*Qplus*fplus)[0]+(fminus.T*Qminus*fminus)[0]
zero('SF11 entire parameterized coefficient identity',gram-target)
for z in x:
    zero('SF11 full cancellation of auxiliary '+str(z),gram.diff(z))

coordinates=[]
columns=[]
for block, ff in enumerate([fplus,fminus]):
    for i in range(len(ff)):
        for j in range(i,len(ff)):
            coordinates.append((block,i,j))
            columns.append(s.Poly(ff[i]*ff[j]*(1 if i==j else 2),u,S,X,Y))
mons=sorted(set().union(*(p.as_dict() for p in columns)))
G=s.Matrix([[p.coeff_monomial(mon) for p in columns] for mon in mons])
require('SF11 exact unrestricted Gram-coordinate count',G.cols==25)
require('SF11 exact independent coefficient count',G.rows==16)
require('SF11 exact coefficient-map rank',G.rank()==16)
qvec=s.Matrix([(Qplus if k==0 else Qminus)[i,j] for k,i,j in coordinates])
Aux=qvec.jacobian(x)
require('SF11 nine actual independent auxiliary directions',Aux.rank()==9)
matrix_zero('SF11 actual auxiliaries span the entire homogeneous kernel',G*Aux)

# Literal zero-strength faces have exact positive Gram witnesses.
# The four coordinates here are n,m,C,A; their repeated physical events
# are precisely those of the reciprocity-fixed sector.
event_witnesses={
    'n':([0,1,0,0,0],[0,0,0,2,0,0,0,0,0],
         s.diag(0,2,0,0),s.diag(1,0,0,1,0)),
    'm':([0,0,1,0,0],[-2,2,-2,2,-2,2,2,-2,0],
         2*s.Matrix([1,-1,1,-1])*s.Matrix([1,-1,1,-1]).T,s.zeros(5)),
    'C':([0,0,0,1,0],[0,0,0,0,0,0,2,-2,0],
         2*s.Matrix([0,0,-1,1])*s.Matrix([0,0,-1,1]).T,
         2*s.Matrix([0,-1,1,0,0])*s.Matrix([0,-1,1,0,0]).T),
    'A':([0,0,0,0,1],[0,0,0,0,0,0,2,0,0],
         s.diag(0,0,2,0),s.diag(0,1,0,0,1)),
}
for name,(vals,xx,qp,qm) in event_witnesses.items():
    sub={**dict(zip([d,n,m,C,A],vals)),**dict(zip(x,xx))}
    matrix_zero('literal zero-strength event '+name+' Qplus witness',Qplus.subs(sub)-qp)
    matrix_zero('literal zero-strength event '+name+' Qminus witness',Qminus.subs(sub)-qm)

# Export the evaluated pencil in the same upper-triangular COO convention
# as the archived C2 dictionaries.  The scalar blocks enforce all orthants
# and both exact equality constraints in original seven coordinates.
variables=list(source)+list(x)
block_matrices=[Qplus,Qminus]+[s.Matrix([[z]]) for z in
    [d,n,m,C,A,D-m,m-D,alpha-C,C-alpha]]
blocks=[]
entries=[]
offset=0
for bi,block in enumerate(block_matrices):
    blocks.append({'index':bi,'offset':offset,'order':block.rows})
    for i in range(block.rows):
        for j in range(i,block.cols):
            val=s.expand(block[i,j])
            constant=val.subs(dict.fromkeys(variables,s.Integer(0)))
            if constant:
                entries.append([0,offset+i,offset+j,str(constant)])
            remainder=val-constant
            for vi,var in enumerate(variables,1):
                fac=val.diff(var)
                require('affine rational COO entry '+str((bi,i,j,vi)),
                        fac.is_Rational is True)
                if fac:
                    entries.append([vi,offset+i,offset+j,str(fac)])
                remainder-=var*fac
            zero('affine entry reconstruction '+str((bi,i,j)),remainder)
    offset+=block.rows
require('SF2 total original-coordinate order 18',offset==18)
require('SF2 seven visible and nine auxiliary coordinates',len(variables)==16)
require('SF2 maximum semidefinite block order 5',max(b['order'] for b in blocks)==5)
require('SF2 exact orthant/equality scalar blocks',
        [b[0,0] for b in block_matrices[2:]]==[d,n,m,C,A,D-m,m-D,alpha-C,C-alpha])
coo={
    'family':'C2 reciprocity-fixed sector',
    'role':'EXACT_AFFINE_LIFT',
    'scope':'Exactly the entire physical C2 cone intersected with D=m and alpha=C; all nonnegative coefficient faces included. No converse for unrestricted C2 is asserted.',
    'variables':list(map(str,variables)),
    'visible_count':7,'auxiliary_count':9,'matrix_order':18,
    'blocks':blocks,'entries':entries,
    'coefficient_convention':'upper triangle; variable index 0 is constant; variables are 1-based in listed order; mirror each off-diagonal once',
    'features':[['S','u*X','u*Y','u^2*S'],['X','Y','u*S','u^2*X','u^2*Y']],
    'state_change':'W=S+T, Z=S-T; X and Y unchanged',
    'testing_reduction':'v=1/u, then all real u by exact sign congruence; SF3-SF7 prove sufficiency on precisely this sector',
    'analytic_input':{'name':'Matrix Fejer-Riesz theorem on the real line with degree bound',
        'source':'Aljaz Zalar, Matrix Fejer-Riesz theorem with gaps, arXiv:1503.06034v1, Theorem 1.2, page 1',
        'url':'https://arxiv.org/abs/1503.06034',
        'application':'A real symmetric degree-four matrix positive semidefinite on the whole real line has real polynomial SOS rows of degree at most two.'},
    'gram_coordinate_count':25,'coefficient_rank':16,'free_dimension':9,
}
export_path=HERE/'C2_RECIPROCAL_FIXED_SECTOR.coo.json'
write_json(export_path,coo)
reloaded=json.loads(export_path.read_text())
reconstructed=s.zeros(18)
for vi,i,j,fac in reloaded['entries']:
    val=s.Rational(fac)*(1 if vi==0 else variables[vi-1])
    reconstructed[i,j]+=val
    if i!=j: reconstructed[j,i]+=val
matrix_zero('reloaded literal COO exactly equals the printed affine pencil',
            reconstructed-s.diag(*block_matrices))

proof=HERE/'c2_reciprocal_fixed_sector.tex'
receipt={
    'status':'PASS_EXACT_FIXED_SECTOR_IDENTITIES_AND_EVALUATED_LIFT',
    'assertion_count':len(CHECKS),'checks':CHECKS,
    'scope':'Complete physical C2 reciprocity-fixed sector D=m, alpha=C. Universal C2 and enlarged K2 coverage remain unproved.',
    'proof':{'path':proof.name,'sha256':sha(proof)},
    'checker':{'path':Path(__file__).name,'sha256':sha(Path(__file__))},
    'export':{'path':export_path.name,'sha256':sha(export_path)},
    'dimensions':{'visible_coordinates':7,'independent_sector_coordinates':5,
                  'auxiliary_variables':9,'largest_block':5,
                  'matrix_order_original_coordinates':18,
                  'matrix_order_independent_coordinates':14},
    'analytic_external_input':coo['analytic_input'],
    'excluded_claims':['Universal C2 coverage','Old or enlarged two-angle Gram exactness','Minimality of the lift size'],
}
write_json(HERE/'check_c2_reciprocal_fixed_sector.json',receipt)
print(json.dumps({k:receipt[k] for k in ['status','assertion_count','dimensions']},indent=2))
