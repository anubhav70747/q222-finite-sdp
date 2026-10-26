#!/usr/bin/env python3
"""Exact sparse C2 corner-pairing audit; no solver, sample search, or network.

Input: saved exact support-face 154-feature map and 24-target row labels.
Claims: (1) 14 integer-corner contact roots give a 78x91 integer map
of rank72; off-target constraints rank54 and the on-target image rank18;
(2) 19 individually support-preserving pair differences form three K4
components and an isolated edge, rank17;
(3) the exact one-dimensional missing coefficient is
[u^2 Y^2] - [v^2 X^2] = 2c^4(b+r)s(s+a^2 c^2t)(h-a^2r)/Delta,
using the complete native C2 parent and exact native fractions.

This is an AFFINE result, not a PSD cone admission! The rank72
14-feature map has no claimed PSD right inverse, and rank17 applies
only to the indicated individually-support-preserving edge squares.
"""
import json, pathlib, itertools
import sympy as sp
import c2_fixed_map as q
HERE=pathlib.Path(__file__).resolve().parent
m=json.loads((HERE/'C2_154_CONTACT_MAP.json').read_text())
target=json.loads((HERE/'C2_CLEARED_POLYNOMIAL_TARGET.json').read_text())
rows=list(map(tuple,m['coefficient_rows']));rid={v:i for i,v in enumerate(rows)}
active=sorted(map(int,target['nonzero_polynomial_coefficients_by_row']))
actset=set(active);active_idx={v:i for i,v in enumerate(active)}
bl=m['support_blocks'];B0=bl[0]
verts=[]
for v in B0:
 i,j,l=v
 if i in (0,6) and j in (0,6) and rid.get((2*i,2*j,l,l)) in actset:
  verts.append(tuple(v))
assert len(verts)==14

def prod(x,y):
 return rid[(x[0]+y[0],x[1]+y[1],min(x[2],y[2]),max(x[2],y[2]))]

# Column convention: symmetric entries are upper triangular; off-diagonal
# coefficients receive the factor 2. Each root is z_i-z_p.
p=verts[0]; vv=verts[1:]; cols=[];names=[]
for j,zj in enumerate(vv):
 for i in range(j+1):
  zi=vv[i];d={}
  for x,y,sgn in ((zi,zj,1),(zi,p,-1),(p,zj,-1),(p,p,1)):
   ix=prod(x,y);d[ix]=d.get(ix,0)+sgn*(1 if i==j else 2)
  cols.append({k:v for k,v in d.items() if v});names.append((i,j))
A=sp.zeros(len(rows),len(cols))
for j,col in enumerate(cols):
 for i,v in col.items(): A[i,j]=v
nonzero=sorted({i for col in cols for i in col});outside=sorted(set(nonzero)-actset)
assert len(nonzero)==78 and len(outside)==54 and len(cols)==91
L=A.extract(outside,list(range(A.cols)))
ker=L.nullspace();K=sp.Matrix.hstack(*ker)
assert K.shape==(91,37)
F=A.extract(active,list(range(A.cols)))*K
assert L.rank()==54 and A.extract(nonzero,list(range(A.cols))).rank()==72
assert F.rank()==18
N=sp.Matrix(m['contact_annihilators'])[:,active]
assert N.rank()==6 and N*F==sp.zeros(6,F.cols)
print('PASS 14-corner 78x91 integer map rank 72, off-target rank54; supported image rank18',flush=True)

# Individually supported pair-square graph. No hidden off-target cancellations.
edges=[];edgecols=[]
for ii,aa in enumerate(verts):
 for jj in range(ii+1,len(verts)):
  bb=verts[jj]
  if prod(aa,bb) in actset:
   e=sp.zeros(len(active),1)
   e[active_idx[prod(aa,aa)]]+=1
   e[active_idx[prod(bb,bb)]]+=1
   e[active_idx[prod(aa,bb)]]-=2
   edges.append((ii,jj));edgecols.append(e)
assert len(edges)==19
E=sp.Matrix.hstack(*edgecols)
assert E.rank()==17 and N*E==sp.zeros(6,E.cols)
# Each edge is part of one of three disjoint K4s or the isolated K2.
components=[];remaining=set(range(len(verts)))
while remaining:
 frontier={min(remaining)}; seen=set()
 while frontier:
  x=frontier.pop();seen.add(x)
  frontier.update(y for y in remaining-seen if (min(x,y),max(x,y)) in edges)
 components.append(sorted(seen));remaining-=seen
assert sorted(map(len,components))==[2,4,4,4]
print('PASS 19 exactly supported contact-square edges: 3 K4 + 1 K2, rank17',flush=True)

# Difference of two active target rows annihilates all 19 simple edges.
ky=(12,0,2,2);kx=(0,12,1,1)
extra=sp.zeros(1,len(active))
extra[0,active_idx[rid[ky]]]=1;extra[0,active_idx[rid[kx]]]=-1
assert extra*E==sp.zeros(1,E.cols) and sp.Matrix.vstack(N,extra).rank()==7
assert (extra*F)!=sp.zeros(1,F.cols)
print('PASS only missing simple-edge direction: coefficient [U^2 phi_Y^2]-[V^2 phi_X^2]',flush=True)

# Identity for every native point in Omega. Source weights symbolic first,
# then exact native fractions. Source contact scalar multipliers are positive.
psi=(q.b/q.a,-q.c*q.r,q.c*(q.b/q.a-q.a*q.t),sp.Integer(1))
pco=q.target_template()
def norm(key):
 ix,jx,i,j=key
 return sp.factor(pco[key]*q.c**(ix//6)*(q.a*q.c)**(jx//6)*psi[i]*psi[j])
left=sp.factor(norm(ky)-norm(kx))
assert sp.simplify(left-(4*q.c**4*q.A*(q.b-q.a**2*q.t)**2/q.a**2-2*q.a**2*q.c**4*q.r**2*(2*q.n+1)))==0
wt,g=q.native_weights();s=1+q.c**2*q.r;h=q.b-q.a**2*q.t
right=2*q.c**4*(q.b+q.r)*s*(s+q.a**2*q.c**2*q.t)*(h-q.a**2*q.r)/g['Delta']
assert sp.cancel((left-right).subs(wt,simultaneous=True))==0
assert sp.cancel((2*q.A*h-q.a**2*q.r*(2*q.n+1)).subs(wt,simultaneous=True))==0
assert sp.cancel((left-2*q.c**4*q.r*(2*q.n+1)*(h-q.a**2*q.r)).subs(wt,simultaneous=True))==0
print('PASS exact native gap factor: difference = positive prefactor*(h-a^2*r)',flush=True)

# Deterministic integer right inverse on the 18-dimensional supported
# contact-stationary coefficient space. This does not assert PSD.
columns=list(F.rref()[1]);assert len(columns)==18
BF=F[:,columns]
independent_rows=list(BF.T.rref()[1]);assert len(independent_rows)==18
Q=BF.extract(independent_rows,range(18))
assert Q.det()!=0
R=K[:,columns]*Q.inv()
V=F[:,columns]*Q.inv()
assert V.extract(independent_rows,range(18))==sp.eye(18)
assert L*R==sp.zeros(L.rows,R.cols)
assert sp.Matrix(m['contact_annihilators'])[:,active]*F[:,columns]==sp.zeros(6,18)
# The 24x18 map V reproduces all 24 rows whenever the six contact
# conditions hold. Check that it is identity on the full 18-dimensional
# target subspace ker N, not only on isolated sources.
for col in F.columnspace():
 assert V*col[independent_rows,:]==col
print('PASS deterministic 18x18 rational RIGHT INVERSE for all native stationary targets, AFFINE ONLY',flush=True)
if __name__=='__main__':
 # Save only original exact rational matrix/data; no runtime pickle.
 def sparse(v):
  return [[i,str(val)] for i,val in enumerate(v) if val]
 obj={'schema':'c2-corner-affine-18-not-psd','14_vertices':verts,
      'integer_91_column_indices':names,'active_24_rows':active,
      'active_rows_keys':[rows[i] for i in active],
      '19_individually_supported_edges':edges,'components':components,
      'offtarget_54_rows':outside,
      'retained_coeff_row_indices':independent_rows,
      'rational_right_inverse_91x18_sparse_columns':[sparse(R[:,j]) for j in range(R.cols)],
      'IMPORTANT':'Affine exact only. Native strict C2 L6 PSD admission is unproved.'}
 (HERE/'C2_CORNER_AFFINE_18.json').write_text(json.dumps(obj,indent=2)+'\n')
 print('WROTE C2_CORNER_AFFINE_18.json',flush=True)
# An exact affine decomposition of *every* source-stationary 24-row target:
# 17 individually supported two-feature contact-square channels plus one
# contact-zero (generally indefinite) 14-feature bridge. Positivity is NOT asserted.
bridge_candidates=[j for j in range(F.cols) if (extra*F[:,j])[0]!=0]
assert bridge_candidates
jbridge=bridge_candidates[0]
tbridge=sp.Rational(1)/(extra*F[:,jbridge])[0]
bridge_vec=tbridge*F[:,jbridge]
bridge_gram=tbridge*K[:,jbridge]
assert (extra*bridge_vec)[0]==1 and L*bridge_gram==sp.zeros(L.rows,1)
ind_edges=list(E.rref()[1]);assert len(ind_edges)==17
Ebas=E[:,ind_edges]
ind_rows=list(Ebas.T.rref()[1]);assert len(ind_rows)==17
J=Ebas.extract(ind_rows,range(17))
assert J.det()!=0
RE=J.inv()
# Verify each of 18 basis vectors of ker N is reconstructed exactly.
for col in F.columnspace():
 residual=col-(extra*col)[0]*bridge_vec
 coeff=RE*residual[ind_rows,:]
 assert Ebas*coeff==residual
print('PASS 18-dimensional exact affine bridge identity: 17 contact-square edge channels + 1 source-algebraic bridge',flush=True)
if __name__=='__main__':
 def sparse(v):
  return [[i,str(z)] for i,z in enumerate(v) if z]
 obj={
  'schema':'c2-24-target-17-edge-plus-bridge-affine-not-psd',
  'bridge_gram_91_contact_reduced_upper_triangular_sparse':sparse(bridge_gram),
  'bridge_coeff_24':list(map(str,bridge_vec)),
  '17_edge_indices_into_19':ind_edges,
  '17_selected_coefficient_rows_into_24':ind_rows,
  'edge_coefficient_inverse_17x17':[[str(RE[i,j]) for j in range(RE.cols)] for i in range(RE.rows)],
  '19_edges':edges,
  'active_24_row_keys':[rows[i] for i in active],
  'normalized_source_extra_coefficient':'2*c^4*(b+r)*(1+c^2*r)*(1+c^2*r+a^2*c^2*t)*(b-a^2*t-a^2*r)/Delta',
  'PSD_STATUS':'UNPROVED: the bridge is not asserted PSD, edge weights may be signed, and extra fractional blocks remain indispensable for the known L2-excluded C2 source.'
 }
 (HERE/'C2_AFFINE_BRIDGE_17_PLUS_1.json').write_text(json.dumps(obj,indent=2)+'\n')
 print('WROTE C2_AFFINE_BRIDGE_17_PLUS_1.json')
# Exact Schur compression of the bridge on the nullspace of the 19-edge
# contact-square network. This is the first failure of that SPECIFIC
# positive-pairing architecture, not an L6 obstruction.
red=sp.zeros(len(verts)-1)
for entry,value in sparse(bridge_gram):
 i,j=names[entry];red[i,j]=red[j,i]=sp.Rational(value)
Vmat=sp.Matrix.vstack(-sp.ones(1,len(verts)-1),sp.eye(len(verts)-1))
Hbridge=Vmat*red*Vmat.T
Qind=sp.zeros(len(verts),len(components))
for k,comp in enumerate(components):
 for i in comp:Qind[i,k]=1
Hquot=Qind.T*Hbridge*Qind
expected=sp.Matrix([[1,0,-1,0],[0,-1,0,1],[-1,0,1,0],[0,1,0,-1]])
assert Hquot==expected
assert sp.factor(Hquot.charpoly().as_expr())==sp.Symbol('lambda')**2*(sp.Symbol('lambda')-2)*(sp.Symbol('lambda')+2)
for i,j in edges:
 vec=sp.zeros(len(verts),1);vec[i]=1;vec[j]=-1
 assert Qind.T*vec==sp.zeros(len(components),1)
print('PASS bridge projected to 4 component constants: eigenvalues (+2,-2,0,0); all 19 simple edge squares vanish there',flush=True)
# As a function of source, Delta(theta) changes sign on actual strict C2.
# These two previously published rational native sources are used only
# to establish sign change, not for any feasibility claim.
points=[(sp.Rational(13,600),sp.Rational(101,1000),sp.Rational(47,15),sp.Rational(27,35),sp.Rational(13,154)),
        (sp.Rational(2),sp.Rational(2),sp.Rational(1,4),sp.Rational(1,4),sp.Rational(1,2))]
vals=[]
for pt in points:
 assert q.is_strict_native(pt)
 vv0=dict(zip((q.a,q.b,q.c,q.t,q.r),pt))
 vals.append(sp.sign((h-q.a**2*q.r).subs(vv0)))
assert vals==[1,-1]
print('PASS h-a^2*r has both signs on the genuine strict native domain Omega',flush=True)
