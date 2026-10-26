"""Exact independent rational/algebraic audit of EN10 reflection and fixed feature transport.
No numerical experiments or optimization.
"""
import sympy as sp

# Distinct symbols for Bell coefficients, testing angles, unrestricted real state.
d, lam, n, m, C, A, D, alpha, beta=sp.symbols('d lam n m C A D alpha beta')
u,v,W,X,Y,Z=sp.symbols('u v W X Y Z', nonzero=True)
coeffs=(d,lam,n,m,C,A,D,alpha,beta)
ref=(d,n,lam,C,m,alpha,beta,A,D)
phi=(W,X,Y,Z)
G=sp.Matrix([[v,-1,0,0],[-1,-v,0,0],[0,0,-v,1],[0,0,1,v]])
assert G.T*G==(1+v**2)*sp.eye(4)

def parent(c,angle_u,angle_v,state):
    d,lam,n,m,C,A,D,alpha,beta=c
    u,v=angle_u,angle_v
    W,X,Y,Z=state
    Du,Dv=1+u**2,1+v**2
    return (4*lam*Du*(W+v*X)**2 + 4*n*Du*Dv*X**2
      +4*m*(v*W-X+u*v*Y-u*Z)**2 +4*C*Dv*(u*W-Y)**2
      +4*A*Du*Dv*Y**2 +(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
      +(4*alpha+2*d)*Du*(v*Y-Z)**2
      +(4*beta+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z**2)

F=sp.expand(parent(coeffs,u,v,phi))
Fref=parent(ref,1/u,v,G*sp.Matrix(phi))
assert sp.cancel(sp.expand(u**2*Fref-(1+v**2)*F))==0
print('PASS exact SG2 full 4-state Bell/source permutation and EN10 congruence')

# Exact original native C1 source transformation, using positive real contact variables.
k,f,g,a,t=sp.symbols('k f g a t',nonzero=True)
kappa=k*k; eta=f*f;zeta=g*g; r=1-a*eta
th=t*zeta; h=1-th
H=1+(1-a)*eta
E=1-kappa*eta*h
Q=1-kappa*zeta*r
T0=1+kappa*h
R0=1+kappa*r
M=R0-t
J=T0-a
G0=h+zeta
V0=kappa*eta*zeta

uc=k; vc=k*f*g; psi=sp.Matrix([g/f,-k*h,k*g*r/f,1])
psi_ref=G.subs(v,vc)*psi
zprime=sp.factor(psi_ref[3]); normalized=sp.simplify(psi_ref/zprime)
# EN13 native inverse gives k',f'^2,g'^2,a',t' (no square-root choices).
kp=1/k
etap=sp.factor(vc/(kp*normalized[0]));zetap=sp.factor(normalized[0]*vc/kp)
rp=sp.factor(normalized[2]/(kp*normalized[0]));thetap=sp.factor(1+normalized[1]/kp)
ap=sp.factor((1-rp)/etap);tp=sp.factor(thetap/zetap)
expected=(kappa*zeta*H/G0,kappa*eta*G0/H,Q/G0,eta*J/H,M/(kappa*H),J/(kappa*G0))
for given,goal in zip((etap,zetap,rp,thetap,ap,tp),expected):
    assert sp.cancel(given-goal)==0,sp.factor(given-goal)
print('PASS transformed native kappa eta zeta r theta a t and state normalizer')

kp2=1/kappa
rp0=Q/G0;hp=E/H
eta_p=kappa*zeta*H/G0; zeta_p=kappa*eta*G0/H
ap0=M/(kappa*H);tp0=J/(kappa*G0)
trans = {
'Hnew':(rp0+eta_p,(1+V0)/G0),
'Gnew':(hp+zeta_p,(1+V0)/H),
'Enew':(1-kp2*eta_p*hp,h*(1+V0)/G0),
'Qnew':(1-kp2*zeta_p*rp0,r*(1+V0)/H),
'Mnew':(1+kp2*rp0-tp0,a*(1+V0)/(kappa*G0)),
'Jnew':(1+kp2*hp-ap0,t*(1+V0)/(kappa*H)),
'c3new_gap':(1-kp2*eta_p+tp0*(kp2*eta_p*zeta_p),h*(1+V0)/G0),
'Vnew':(kp2*eta_p*zeta_p,V0),
}
for label,(left,right) in trans.items():
    assert sp.cancel(left-right)==0,(label,sp.factor(left-right))
print('PASS all eight reflected native positive-gap identities')
assert sp.cancel(zprime-k*g*H/f)==0
print('PASS positive source-only state contact normalizer Zprime = k*g*H/f')

# Z16 fixed reflected features.  Use formal U,V and ordinary source state placeholders.
U,V=sp.symbols('U V',positive=True)
WW,XX,YY,ZZ=sp.symbols('WW XX YY ZZ')
S=sp.Matrix([WW,XX,YY,ZZ]); Sr=G.subs(v,V)*S
z16=lambda p,q,ss:sp.Matrix([p*ss[0],q*ss[0],ss[1],p*q*ss[1],ss[2],p*q*ss[2],p*ss[3],q*ss[3],p*q*ss[0],q*ss[1],p*ss[2],q*ss[2],ss[3],p*q*ss[3],sp.sqrt(p*q)*ss[0],sp.sqrt(p*q)*ss[3]])
zr=sp.simplify(U*z16(1/U,V,Sr))
assert len(zr)==16
# Powers of U and V in each individual row: 0, 1/2 or 1 in U; max 2 in V.
for expr in zr:
    ex=sp.expand(expr)
    for term in sp.Add.make_args(ex):
        term=sp.powsimp(term,force=True)
        uu=term.as_powers_dict().get(U,0); vv=term.as_powers_dict().get(V,0)
        assert 0<=uu<=1 and 0<=vv<=2 and int(6*uu)==6*uu and int(6*vv)==6*vv,(expr,term,uu,vv)
print('PASS all 16 fixed reflected algebraic rows lie in the 7x13 sixth-root grid')

# Double orientation is an involution as coefficient permutation and projective angle transform.
assert (d,lam,n,m,C,A,D,alpha,beta)==(ref[0],ref[2],ref[1],ref[4],ref[3],ref[7],ref[8],ref[5],ref[6])
print('PASS source reflection is an involutive integral permutation')

if __name__=='__main__':
    print('16 explicit unmultiplied reflected rows:')
    for i,row in enumerate(zr,1):print(f'{i:2}  {sp.sstr(row)}')

# The fixed SG1 linear map from canonical source coefficients to conventional
# affine Bell coefficients [1,A0,A1,B0,B1,E00,E01,E10,E11].
E_Bell=sp.zeros(9,9)
E_Bell[0,0]=1
for row,scalar in [(5,-sp.Rational(1,2)),(6,sp.Rational(1,2)),(7,-sp.Rational(1,2)),(8,-sp.Rational(1,2))]:E_Bell[row,0]=scalar
for col,(x,y,sign_a,sign_b) in enumerate([(0,1,1,1),(0,0,1,-1),(1,1,1,-1),(1,0,-1,1),(0,0,-1,1),(1,1,-1,1),(0,1,-1,-1),(1,0,1,-1)],1):
    E_Bell[0,col]=1
    E_Bell[1+x,col]=sign_a
    E_Bell[3+y,col]=sign_b
    E_Bell[5+2*x+y,col]=sign_a*sign_b
Sig=sp.zeros(9,9)
for row,col in enumerate([0,2,1,4,3,7,8,5,6]):Sig[row,col]=1
assert E_Bell.det()!=0
B_std=sp.simplify(E_Bell*Sig*E_Bell.inv())
assert B_std*B_std==sp.eye(9)
assert B_std==sp.Matrix([[1,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0],[0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1,0]])
print('PASS SG1 source basis-to-standard Bell linear map, 9x9, and involution')

# Build the fixed 364 x 16 signed sparse feature transport, containing only 0,+1,-1.
labels=('W','X','Y','Z')
base16=[(6,0,0),(0,6,0),(0,0,1),(6,6,1),(0,0,2),(6,6,2),(6,0,3),(0,6,3),(6,6,0),(0,6,1),(6,0,2),(0,6,2),(0,0,3),(6,6,3),(3,3,0),(3,3,3)]
ext=[(i,j,l) for i in range(7) for j in range(13) for l in range(4)]
lookup={w:i for i,w in enumerate(ext)}
G_terms={0:((6,0,1),(0,1,-1)),1:((0,0,-1),(6,1,-1)),2:((6,2,-1),(0,3,1)),3:((0,2,1),(6,3,1))}
L_sparse=[]
for (i,j,l) in base16:
    block=[]
    for shift,ell,sign in G_terms[l]:
        block.append((lookup[(6-i,j+shift,ell)],sign))
    L_sparse.append(block)
assert len(ext)==364 and len(L_sparse)==16 and sum(map(len,L_sparse))==32
# Check expansion reproduces the actual exact h16 row expressions.
state_syms=(WW,XX,YY,ZZ)
for h_expr,block in zip(zr,L_sparse):
    ans=sum(sgn*U**sp.Rational(ext[idx][0],6)*V**sp.Rational(ext[idx][1],6)*state_syms[ext[idx][2]] for idx,sgn in block)
    assert sp.simplify(h_expr-ans)==0
print('PASS 32 signed sparse entries reproduce all 16 reflected features on the fixed 364 grid')

# The two canonical embeddings J0, J1 of 196 rows into the 364 grid cover every
# extended coordinate, so Q_ext=J0*Q*J0.T + J1*Q*J1.T is PD whenever Q is PD.
base196=[(i,j,l) for i in range(7) for j in range(7) for l in range(4)]
j0=[lookup[(i,j,l)] for i,j,l in base196]
j1=[lookup[(i,j+6,l)] for i,j,l in base196]
assert set(j0)|set(j1)==set(range(364)) and len(j0)==len(j1)==196
print('PASS full 196-row to 364-row two-embedding coverage for Dv multiplication and strict Slater')

# Native coefficient formula from SD1/EN13, evaluated exactly (no numerical sweep).
def native_c1(k2, et, ze, a, t):
    Vc=k2*et*ze; Lc=1-a*et; Hc=1+(1-a)*et; Qc=1-k2*ze+a*Vc
    R=1+k2*Lc; T=1+k2-k2*ze*t; M=R-t; J=T-a
    P4=(1+k2+k2*et*(a-k2-1)-ze*k2*t
        +Vc*(-ze*k2*t**2+2*k2*t+a*(1-t)+t))
    Pal=(1+k2-a*et*k2+ze*((1+k2)-(2+k2)*t)
         +Vc*(-ze*t**2+a*(t-1)+t))
    Pd=((1+k2)**2-a*et*k2**2 +ze*k2*((k2+1)-(k2+2)*t)
        +Vc*(-ze*k2*t**2+k2*a*(t-1)+k2*t+k2+1))
    c2=2*(1+k2)*(1-t*ze)*Lc*(1+Vc)*Hc*Qc
    cb=2*(1+k2)*Lc*(1+(1-t)*ze)*Hc*Qc
    c3=2*(1+k2)*Lc*Hc*Qc*(1-k2*et+t*Vc)
    c5=2*(1+k2)*Lc*Hc*Qc*(1+Vc)
    c4=2*Qc*P4
    ca=2*Lc*Hc*Pal
    cd=4*Lc*Qc*Pd
    dB=2*a*et*cd/c3
    mB=2*a*(k2*et/ze)*c2/(M*c3)
    AB=2*a*et**2*c4/c3
    alphaB=2*a*(ze/(k2*et))*k2**3*et**2*ca/c3
    DB=2*a*c5/(J*c3)-dB/2
    betaB=2*a*(k2*et/ze)*cb/(t*c3)-dB/2
    return tuple(map(sp.factor,(dB,sp.Integer(0),sp.Integer(0),mB,sp.Integer(2),AB,DB,alphaB,betaB)))

# The published exact strict k=4 obstruction witness is useful *only* as a
# source-map consistency check: it is not used as a PSD admission argument.
sub_w={k:sp.Integer(4),f:sp.Rational(1,1000),g:sp.Rational(1,8),a:sp.Integer(1),t:sp.Integer(1)}
c_w=native_c1(sp.Integer(16),sp.Rational(1,1000000),sp.Rational(1,64),sp.Integer(1),sp.Integer(1))
kp2_w=sp.factor((1/kappa).subs(sub_w))
eta_p_w=sp.factor((kappa*zeta*H/G0).subs(sub_w))
zeta_p_w=sp.factor((kappa*eta*G0/H).subs(sub_w))
a_p_w=sp.factor((M/(kappa*H)).subs(sub_w))
t_p_w=sp.factor((J/(kappa*G0)).subs(sub_w))
c_or=native_c1(kp2_w,eta_p_w,zeta_p_w,a_p_w,t_p_w)
# Native source always fixes C=2; the reflected Bell source has C'=m_old.
scale=sp.cancel(2/c_w[3]); reflected=tuple(c_w[j] for j in [0,2,1,4,3,7,8,5,6])
for i,(left,right) in enumerate(zip(c_or,reflected)):
    assert sp.cancel(left-scale*right)==0,(i,left,scale*right)
print('PASS all nine exact native Bell weights under reflection, with essential C=2 renormalization')

# Reproducible, integer sparse transports: writing this file is deterministic.
import json
from pathlib import Path
mapping={
    'source_order':['d','lambda','n','m','C','A','D','alpha','beta'],
    'source_reflection_output_indices':[0,2,1,4,3,7,8,5,6],
    'bell_affine_order':['1','A0','A1','B0','B1','E00','E01','E10','E11'],
    'sg1_to_affine':[[str(x) for x in row] for row in E_Bell.tolist()],
    'affine_reflection':[[str(x) for x in row] for row in B_std.tolist()],
    'state_G_v':[['v','-1','0','0'],['-1','-v','0','0'],['0','0','-v','1'],['0','0','1','v']],
    'monomial_index_convention':'(i,j,l) -> u**(i/6)*v**(j/6)*state[l], zero-based l in W,X,Y,Z',
    'full_extended_grid_shape':[7,13,4],
    'extended_features':[[i,j,l] for i,j,l in ext],
    'oriented_Z16_features':[[i,j,l] for i,j,l in base16],
    'reflected_16_sparse_columns':[[[idx,sign] for idx,sign in block] for block in L_sparse],
    'full_Z6_shape':[7,7,4],
    'Z6_to_ext_unmultiplied_indices':j0,
    'Z6_to_ext_times_v_indices':j1,
    'coefficient_key_shape':[13,25,10],
    'coefficient_linear_map':'coefficient identity (1+v**2)*SG2(c)=Zext.T * Q * Zext; include offdiag Gram factor 2',
    'native_contact_reflection':{
        'kprime':'1/k','eta_prime':'kappa*zeta*H/G','zeta_prime':'kappa*eta*G/H',
        'aprime':'M/(kappa*H)','tprime':'J/(kappa*G)',
        'normalized_contact_state':['f*G/(g*H)','-E/(k*H)','f*Q/(k*g*H)','1'],
        'contact_state_divisor':'k*g*H/f',
        'C2_normalization_of_reflected_C1_source':'(2/m)*sigma(c) if original C=2',
        'definitions':'r=1-a*eta, h=1-t*zeta, H=r+eta, G=h+zeta, E=1-kappa*eta*h, Q=1-kappa*zeta*r, M=1+kappa*r-t, J=1+kappa*h-a'
    },
    'checked_assertion_count_groups':8
}
json_path=Path(__file__).with_name('transport_sparse_maps.json')
json_path.write_text(json.dumps(mapping,indent=2,sort_keys=True)+'\n',encoding='utf8')
assert len({tuple(z) for z in mapping['extended_features']})==364
print('PASS deterministic integer sparse map JSON regenerated')
