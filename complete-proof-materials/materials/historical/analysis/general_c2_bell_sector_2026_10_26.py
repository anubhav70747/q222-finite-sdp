"""Exact full-parent Bell blocks and factored native contact check; no solver."""
import json
from pathlib import Path
import sympy as s

u,v,d=s.symbols('u v d',real=True)
Du,Dv=1+u*u,1+v*v
weights=s.symbols('lambda n m C A D alpha beta',real=True)
roots=[s.Matrix(x) for x in [
 [1,v,0,0],[0,1,0,0],[v,-1,u*v,-u],[u,0,-1,0],
 [0,0,1,0],[u,u*v,-1,-v],[0,0,v,-1],[0,1,0,u]]]
factors=[Du,Du*Dv,1,Dv,Du*Dv,1,Du,Dv]
plus=[s.Matrix(x) for x in [
 [1,v],[0,1],[v-u,-1-u*v],[u,1],[0,-1],
 [u-v,1+u*v],[-1,-v],[u,1]]]
minus=[s.Matrix(x) for x in [
 [1,v],[0,1],[u+v,u*v-1],[u,-1],[0,1],
 [u+v,u*v-1],[1,v],[-u,1]]]
# (W,X,Y,Z)^T = B (p1,p2,q1,q2)^T.
B=s.Matrix([[1,0,1,0],[0,1,0,1],[0,-1,0,1],[1,0,-1,0]])/2
for L,lp,lm in zip(roots,plus,minus):
    assert (2*B.T*L-s.Matrix.vstack(lp,lm)).applyfunc(s.expand)==s.zeros(4,1)
eZ=s.Matrix([0,0,0,1])
core=2*(roots[5]*roots[5].T+Du*roots[6]*roots[6].T
        +Dv*roots[7]*roots[7].T-Du*Dv*eZ*eZ.T)
Kp0=s.Matrix([[u*(u-v),u*(1+u*v)],
              [u*(1+u*v),1+v*v+u*u*v*v+u*v]])
Km0=s.Matrix([[u*(u+v),u*(u*v-1)],
              [u*(u*v-1),1+v*v+u*u*v*v-u*v]])
assert (B.T*core*B-s.diag(Kp0,Km0)).applyfunc(s.expand)==s.zeros(4)
Kp,Km,E=d*Kp0,d*Km0,s.zeros(2)
parent=d*core
for w,f,L,lp,lm in zip(weights,factors,roots,plus,minus):
    parent+=4*w*f*L*L.T
    Kp+=w*f*lp*lp.T
    Km+=w*f*lm*lm.T
    E+=w*f*lp*lm.T
assert (B.T*parent*B-s.BlockMatrix([[Kp,E],[E.T,Km]]).as_explicit()).applyfunc(s.expand)==s.zeros(4)
assert s.expand(Km0.det()-u*v*Du*Dv)==0
assert Km0[0,0]==u*(u+v)

# The Schur polynomial identity is checked generically, without expanding
# products of native rational functions or introducing any mixed inverse.
k00,k01,k11=s.symbols('k00 k01 k11')
K=s.Matrix([[k00,k01],[k01,k11]])
ep=s.Matrix(2,2,s.symbols('e0:4'))
P=s.Matrix(2,2,s.symbols('p0:4'))
P=(P+P.T)/2
pvec=s.Matrix(s.symbols('p1 p2')); qvec=s.Matrix(s.symbols('q1 q2'))
F=(pvec.T*P*pvec+2*pvec.T*ep*qvec+qvec.T*K*qvec)[0]
N=K.det()*P-ep*K.adjugate()*ep.T
bvec=K*qvec+ep.T*pvec
assert s.expand(K.det()*F-(bvec.T*K.adjugate()*bvec)[0]-(pvec.T*N*pvec)[0])==0

# Native source comparison and kernel: clear only the common denominator
# after cancelling each evaluated event root against its own source factor.
a,b,c,r,t=s.symbols('a b c r t',positive=True)
z=c*c; Ac=1+z*r; h=b-a*a*t
g=(1+z)*b-(1-r+a*a*z*t); j=1+z*r-t; gap=1-z*h
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b
       +z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
Ns=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
native=[0,((b+r)*Ac*(Ac+a*a*z*t)-r*Delta)/(2*r*Delta),
 (1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Delta),
 Ac*Ns/(2*a*a*t*Delta),a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Delta),
 (1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Delta),
 z*z*(b+r)*Ns/(2*gap*Delta),0]
contact={u:c,v:a*c}
psi=s.Matrix([b/a,-c*r,c*h/a,1])
rv=[(b-a*a*z*r)/a,-c*r,c*g,a*c*t,c*h/a,-a*c*j,-gap,c*(1-r)]
for L,value in zip(roots,rv):
    assert s.expand((L.subs(contact).T*psi)[0]-value)==0
cleared=a*Delta*core.subs(contact)*psi
for w,f,L,value in zip(native,factors,roots,rv):
    if w!=0:
        scalar=s.cancel(a*Delta*4*w*s.sympify(f).subs(contact)*value)
        cleared+=scalar*L.subs(contact)
assert cleared.applyfunc(s.expand)==s.zeros(4,1)
assert s.expand(Ns-((1+z)*(b+r)*h+a*a*(g+j+r*t*(1+a*a*z))))==0
coords=2*B.T*psi
pc=s.Matrix([(b+a)/a,-c*(a*r+h)/a])
qc=s.Matrix([(b-a)/a,c*(h-a*r)/a])
assert (coords-s.Matrix.vstack(pc,qc)).applyfunc(s.expand)==s.zeros(4,1)

# The discriminant is for a*Psi, not the normalized state matrix.
M=s.Matrix([[b,-a*c*r],[c*h,a]])
Gram=M.T*M
disc=(b*b+z*h*h-a*a*(1+z*r*r))**2+4*a*a*z*(h-b*r)**2
assert s.expand(s.trace(Gram)**2-4*Gram.det()-disc)==0
assert s.expand(disc.subs({b:a,t:(1-r)/a}))==0

data={
 'full_canonical_parent_Bell_block_identity_exact':True,
 'core_cross_block_zero':True,
 'core_plus':[[str(w) for w in row] for row in Kp0.tolist()],
 'core_minus':[[str(w) for w in row] for row in Km0.tolist()],
 'event_order':list(map(str,weights)),
 'event_factors':list(map(str,factors)),
 'event_plus_vectors':[[str(w) for w in lp] for lp in plus],
 'event_minus_vectors':[[str(w) for w in lm] for lm in minus],
 'minus_core_determinant':'u*v*(1+u^2)*(1+v^2)',
 'minus_core_first_pivot':'u*(u+v)',
 'positive_block_domain':'d>0,u>0,v>0, all event weights nonnegative',
 'generic_fraction_free_Schur_identity_exact':True,
 'all_four_native_contact_kernel_entries_exact':True,
 'native_event_roots_cancelled_before_common_denominator_expansion':True,
 'contact_Bell_coordinate_formula_exact':True,
 'Schmidt_discriminant_identity_exact':True,
 'Schmidt_discriminant_scale':'formula belongs to a*Psi; divide by a^4 for the unscaled contact matrix',
 'maximal_entanglement_condition':'b=a and a*t=1-r; equivalent to qc=0 in this native sign chart',
 'native_N_positive_identity_exact':True,
 'scope':'Existing pointwise Bell-sector/Schur reduction; no mixed inverse or nonlinear numerator is asserted to define a finite affine SDP.',
 'no_numerical_optimization_or_parameter_sweep':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,indent=2))
