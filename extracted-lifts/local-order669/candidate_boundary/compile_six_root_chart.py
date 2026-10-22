"""Exact fixed post-TF1 chart derived from the current full C1 identity.

This is a LOCAL source-guarded construction, not a full Q222 lift.
Requires SymPy. No numerical SDP or floating-point PSD test is used.
"""
from pathlib import Path
import json
import sys
from fractions import Fraction
from math import isqrt

import sympy as S

HERE = Path(__file__).resolve().parent
EPS = S.Rational(1,4096)
ETA = S.Rational(1,56)
PHI = ['W','X','Y','Z','u*Z','w*Z','u*W','u*w*W','u**2*W','u**3*w*W',
       'w*X','u*w*Z','u*X','u**2*w*W','u*Y','w*Y','u*w*Y']
ALPHAS = [(0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(2,1),(3,0)]
EVENT_GROUPS = [([10,11],[1]),([12,15,7],[2,7]),([4,16,13],[3,8]),([14,11,8],[6,9])]

def source(x,sigma,b,z,h):
    r=b*z; Q=1-b; theta=1-h; e=x*h/sigma; E=1-e; K=1+x
    Rs=1+x*r; M=Rs-x*z*theta; J=1+x*h-sigma*(1-r)
    P=sigma+K+sigma*z*(1+b*x+2*x*h)+sigma*b*x*x*h*z*z+x*x*h*theta*z
    p=(z*(1+x*h)*E+(1-r)*(1+x*h*z))/(2*P)
    ell0=(1+sigma*r)*(e*Rs*(1+sigma*z)+M*E)/(2*P)
    kappa=K*(1+sigma*r)*(1+sigma*z)/(2*M*P)
    H=h/sigma
    m=kappa*x*x*z*H
    C=K*z*(1+sigma*r)*E/(2*(1-r)*P)
    A=p/r; alpha=ell0/(z*Q)
    D=(K*(1+sigma*z)*(1+sigma*r)-J*P)/(2*J*P)
    beta=(K*(1+sigma*r)*(1+x*h*z)-theta*P)/(2*theta*P)
    a=p/b; ell=ell0/Q; g=h-Q
    nu=(2*D-2*beta*theta+h)/sigma
    lam=(8*D*h+4*beta*theta**2+2*h*h)/sigma-4*b*p*sigma*z
    L=lam/sigma; s=S.sqrt(1+L*L); q=S.sqrt(z)
    Bx=8*kappa*e
    Ay=8*(p-ell0)+4*z*((1+2*m)*g+2*D*h+2*p*sigma)
    Rx=-4*(2*D*b+(1+2*m)*g-2*beta*theta)
    R0=-4*(2*(C+D)*(r-1)+2*p+r-1)
    chi=8*(p*b+ell0*Q)-4*h+4*e*h*(2*D+1)-8*D*h+4*z*(4*D*b*h+4*p*b*sigma+(1+2*m)*g*g+2*beta*theta**2)
    k1=2*(1+2*m)*g*g+8*D*h*b+8*p*b*sigma+4*beta*theta**2
    dG=4*D/sigma**2
    return locals()

def matrices(d):
    # Coefficient matrices in roots Nx,Ny,R,Gamma,A_theta,Vhat.
    x,sigma,b,z,h,q,s,L,nu,p,Q,m,D,C,kappa,H=dget(d,'x sigma b z h q s L nu p Q m D C kappa H')
    Bx,Ay,Rx,R0,chi,k1,lam,r=dget(d,'Bx Ay Rx R0 chi k1 lam r')
    E={alpha:S.zeros(6) for alpha in ALPHAS}
    def put(alpha,i,j,v):
        E[alpha][i,j]=v; E[alpha][j,i]=v
    G0=S.diag(2,2,b*b+Q*Q,d['dG'],4*kappa*H,1)
    G0[1,2]=G0[2,1]=2*b-1
    G0[3,5]=G0[5,3]=L/(2*s)
    put((0,0),2,2,z*(8*b*D+8*b*p*sigma+4*b*sigma*(2*C+2*D+1)*(r-1)))
    put((0,0),1,2,2*z*(sigma*(2*C+1)*(r-1)+2*D*(1+sigma*(r-1))+2*p*sigma))
    put((0,0),3,5,-s*x/2)
    put((1,0),0,3,-2*nu)
    put((1,0),1,4,Ay/2)
    put((1,0),2,4,(chi+lam+s*s*sigma*x)/2)
    put((0,1),1,3,q*R0/2)
    put((0,1),2,3,b*q*R0/2)
    put((0,1),0,4,q*Bx/2)
    put((1,1),0,1,-2*q*(2*m+2*D+1))
    put((1,1),0,2,q*Rx/2)
    put((2,0),1,1,2*z*(1+2*m))
    put((2,0),2,2,z*k1)
    put((2,0),1,2,2*z*(2*D*h+2*p*sigma+(1+2*m)*d['g']))
    put((2,0),3,3,L)
    put((2,0),4,4,-s*s*sigma)
    put((2,0),3,5,s/2)
    put((0,2),1,1,2*z*(1+2*D))
    put((0,2),2,2,2*z*b*b*(2*D+1))
    put((0,2),1,2,2*z*b*(2*D+1))
    put((2,1),1,3,-4*p*q)
    put((2,1),2,3,-4*b*p*q)
    put((3,0),2,4,-s*s*sigma/2)
    T=S.zeros(6,10)
    T[0,1]=1; T[0,4]=h
    T[1,2]=1/q; T[1,5]=-b/q
    T[2,5]=1/q; T[2,6]=-q
    T[3,3]=sigma; T[3,7]=-1
    T[4,0]=-q*x; T[4,8]=q
    T[5,7]=-s*x; T[5,9]=s
    return G0,E,T

def dget(d,names): return [d[k] for k in names.split()]

def event_rows(d):
    q,b,Q,x,sigma,h,e=dget(d,'q b Q x sigma h e')
    m,D,C,beta,a,ell,p,z=dget(d,'m D C beta a ell p z')
    def row(items):
        v=S.zeros(17,1)
        for i,c in items: v[i]=c
        return v
    return [
      (4*beta,row([(1,1),(4,h)]),'4 beta Nx^2'),
      (4*beta+2+4*m,row([(10,1),(11,h)]),'(4 beta+2+4m)(w Nx)^2'),
      (2*(2*D+1),row([(12,1),(7,e)]),'2(2D+1)(u X_e)^2'),
      (4*b*p*z,row([(4,sigma),(13,-1)]),'4bpz(u Gamma)^2'),
      (4*a+2*z*(2*C+1)-1,row([(2,1/q),(6,-b*q)]),'clipped J_C excess'),
      (4*ell-1,row([(2,1/q),(5,-1/q),(6,Q*q)]),'clipped J_ell excess'),
      (4*a,row([(14,1/q),(8,-b*q)]),'4a(u J_C)^2'),
      (4*a+4*C*z,row([(15,1/q),(7,-b*q)]),'(4a+4Cz)(w J_C)^2'),
      (4*a,row([(16,1/q),(13,-b*q)]),'4a(uw J_C)^2'),
      (4*ell,row([(14,1/q),(11,-1/q),(8,Q*q)]),'4ell(u J_ell)^2')]

def compile_blocks(d):
    G0,E,T=matrices(d)
    base=G0+E[(0,0)]-S.eye(6)/8
    T17=T.row_join(S.zeros(6,7))
    B0=T17.T*base*T17
    for weight,row,_ in event_rows(d): B0 += weight*row*row.T
    blocks=[B0]
    for alpha in ALPHAS[1:]:
        scaled=EPS**sum(alpha)*E[alpha]
        for sign in [1,-1]: blocks.append(T.T*(ETA*S.eye(6)+sign*scaled)*T/2)
    return [B.applyfunc(S.simplify) for B in blocks]

def compile_structured_blocks(d):
    G0,E,T=matrices(d)
    events=event_rows(d)
    B0=T.T*(G0+E[(0,0)]-S.eye(6)/8)*T
    for j in (0,4,5):
        weight,row,_=events[j]
        assert not any(row[10:,0])
        r=row[:10,0];B0+=weight*r*r.T
    blocks=[B0]
    for alpha in ALPHAS[1:]:
        scaled=EPS**sum(alpha)*E[alpha]
        for sign in [1,-1]:blocks.append(T.T*(ETA*S.eye(6)+sign*scaled)*T/2)
    for indices,events_in_group in EVENT_GROUPS:
        B=S.zeros(len(indices))
        for j in events_in_group:
            weight,row,_=events[j]
            assert all(i in indices for i,x in enumerate(row) if x)
            r=S.Matrix([row[i] for i in indices]);B+=weight*r*r.T
        blocks.append(B)
    return [B.applyfunc(S.simplify) for B in blocks]

def interval_sqrt_rational(q, bits=96):
    q=Fraction(str(q)); scale=1<<bits
    k=isqrt(q.numerator*scale*scale//q.denominator)
    return S.Rational(k,scale),S.Rational(k+1,scale)

def rational_upper_abs(expr, rad):
    # Every occurrence lies in Q(rad), with rad^2 rational. Bounds are exact.
    radicals=[a for a in rad.atoms(S.Pow) if a.exp==S.Rational(1,2)]
    if radicals: rad=radicals[0]
    expr=S.radsimp(expr)
    a=S.simplify((expr+expr.xreplace({rad:-rad}))/2)
    b=S.simplify((expr-expr.xreplace({rad:-rad}))/(2*rad))
    if not (a.is_Rational and b.is_Rational):
        raise ValueError(('not quadratic',expr,a,b))
    lo,hi=interval_sqrt_rational(rad**2)
    low=a+(b*lo if b>=0 else b*hi)
    high=a+(b*hi if b>=0 else b*lo)
    return max(abs(low),abs(high))

def verify_fixture():
    sigma=S.Rational(1,10000)
    d=source(sigma*sigma,sigma,S.Rational(1,2),S.Integer(1),3*sigma/2)
    G0,E,T=matrices(d)
    assert all(d[k]>0 for k in ['E','M','J','m','C','A','alpha','D','beta'])
    assert d['dG']>=1 and 4*d['kappa']*d['H']>=1 and 4*d['a']>=1 and 4*d['ell']>=1
    bounds={}
    for alpha,mat in E.items():
        value=max(sum(rational_upper_abs(mat[i,j],d['s']) for j in range(6)) for i in range(6))*EPS**sum(alpha)
        bound=S.Rational(1,12) if alpha==(0,0) else ETA
        assert value<=bound,(alpha,value,bound)
        bounds[str(alpha)]={'exact_rational_upper':str(value),'required_upper':str(bound)}
    blocks=compile_blocks(d)
    assert len(blocks)==15 and sum(B.rows for B in blocks)==157
    assert all(c.is_Rational for B in blocks for c in B)
    for w,_,_ in event_rows(d): assert w>=0
    u,w,W,X,Y,Z=S.symbols('u w W X Y Z',real=True)
    env={'u':u,'w':w,'W':W,'X':X,'Y':Y,'Z':Z}
    phi=S.Matrix([S.sympify(a,locals=env) for a in PHI])
    actual=(phi.T*blocks[0]*phi)[0]
    for k,alpha in enumerate(ALPHAS[1:]):
        y=(u/EPS)**alpha[0]*(w/EPS)**alpha[1]
        short=phi[:10,0]
        actual += (1+y)*(short.T*blocks[1+2*k]*short)[0]+(1-y)*(short.T*blocks[2+2*k]*short)[0]
    m,C,A,D,ap,beta=dget(d,'m C A D alpha beta')
    Du=1+u*u;Dw=1+w*w
    parent=4*m*(W-w*X+u*Y-u*w*Z)**2+4*C*Dw*(u*W-Y)**2+4*A*Du*Dw*Y*Y+(4*D+2)*(u*w*W+u*X-w*Y-Z)**2+(4*ap+2)*Du*(Y-w*Z)**2+(4*beta+2)*Dw*(X+u*Z)**2-2*Du*Dw*Z*Z
    assert S.Poly(S.expand(actual-parent),u,w,W,X,Y,Z).is_zero
    entries=[]
    for k,B in enumerate(blocks):
        for i in range(B.rows):
            for j in range(i,B.cols):
                if B[i,j]:entries.append([k,i,j,str(B[i,j])])
    data={'status':'Exact source-guarded local full-parent Gram; no global Q222 lift claimed',
          'source_parameters':{k:str(d[k]) for k in ['x','sigma','b','z','h']},
          'source_coefficients':{k:str(d[k]) for k in ['m','C','A','D','alpha','beta']},
          'epsilon':str(EPS),'eta':str(ETA),'base_features':PHI,
          'short_features':PHI[:10],'nonconstant_monomials':[list(a) for a in ALPHAS[1:]],
          'block_sizes':[B.rows for B in blocks],'total_block_order':157,
          'records':'[block, row, column, rational]; upper triangle, mirror once',
          'entries':entries}
    (HERE/'six_root_admitted_fixture.coo.json').write_text(json.dumps(data,indent=2)+'\n')
    receipt={'full_parent_coefficient_identity':True,'all15_blocks_exact_rational':True,
             'PSD_proof':'source guard + symmetric row bound + source congruence; no numerical eigentest',
             'source_guard_bounds':bounds,'source_lambda_sign':int(S.sign(d['lam'])),
             'explicit_blocks':15,'total_order':157,'nonzero_upper_records':len(entries),
             'state_coordinates_all_retained':True,'both_Gamma_normal_couplings_retained':True}
    (HERE/'six_root_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    structured=compile_structured_blocks(d)
    sf=[PHI[:10]]*15+[[PHI[i] for i in indices] for indices,_ in EVENT_GROUPS]
    sactual=0
    for k,mat in enumerate(structured):
        f=S.Matrix([S.sympify(a,locals=env) for a in sf[k]])
        if 1<=k<=14:
            alpha=ALPHAS[1:][(k-1)//2];sgn=1 if k%2 else -1
            weight=1+sgn*(u/EPS)**alpha[0]*(w/EPS)**alpha[1]
        else:weight=1
        sactual+=weight*(f.T*mat*f)[0]
    assert S.Poly(S.expand(sactual-parent),u,w,W,X,Y,Z).is_zero
    sentries=[[k,i,j,str(mat[i,j])] for k,mat in enumerate(structured) for i in range(mat.rows) for j in range(i,mat.cols) if mat[i,j]]
    data2=dict(data)
    data2.update({'block_sizes':[a.rows for a in structured],'total_block_order':161,'entries':sentries,
                  'per_block_unweighted_features':sf,'variant':'15principalblocks oforder10 andfourretainedblocks(2,3,3,3)'})
    (HERE/'six_root_structured_fixture.coo.json').write_text(json.dumps(data2,indent=2)+'\n')
    receipt['structured_variant']={'order':161,'blocks':19,'max_block':10,'records':len(sentries),'full_parent_identity':True}
    (HERE/'six_root_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_guard_bounds'},indent=2))

def export_general_formulas():
    names='x sigma b z h q s L nu p Q m D C beta kappa H Bx Ay Rx R0 chi k1 lam r dG a ell e'
    d={name:S.Symbol(name,real=True) for name in names.split()}
    d['g']=d['h']-d['Q']
    G0,E,T=matrices(d)
    out={'status':'Exact LOCAL source-guarded construction; full-set coverage not asserted',
         'root_order':['N_x','N_y','R','Gamma','A_theta','Vhat'],
         'base_features':PHI,'short_features':PHI[:10],
         'testing_radius':str(EPS),'source_radius_for_lambda_nonpositive':'1/1024',
         'scalar_allocation_per_nonconstant_monomial':str(ETA),
         'source_coefficient_matrix':[[str(x) for x in row] for row in T.tolist()],
         'constant_reference':[[str(x) for x in row] for row in G0.tolist()],
         'coefficient_matrices':{str(alpha):[[str(x) for x in row] for row in mat.tolist()] for alpha,mat in E.items()},
         'retained_event_rows':[{'weight':str(weight),'row':[str(x) for x in row],'identity_term':label} for weight,row,label in event_rows(d)],
         'source_guards':['native P1 domain and nonnegative D,beta','4*a>=1','4*ell>=1','dG>=1','4*kappa*H>=1',
                          'max_i sum_j abs(E_00[i,j]) <= 1/12',
                          'for each nonzero alpha: (1/4096)^sum(alpha)*max_i sum_j abs(E_alpha[i,j]) <= 1/56'],
         'uniform_source_subdomain':'0<x,sigma<=1/1024; 0<z<=1; native P1 gaps; D,beta>=0; lambda<=0',
         'definitions':['q^2=z,q>0','s^2=1+L^2,s>0','L=lambda/sigma','Q=1-b','r=b*z'],
         'block_orders':[17]+[10]*14,
         'normalized_feature_rule':'kappa_t*phi17 for block0; kappa_t*sqrt(1+sign*(u/eps)^a*(w/eps)^b)*phi10 for each alpha and sign',
         'normalizer':'kappa_t>0; kappa_t^2*(1+u^2)*(1+w^2)=1; W^2+X^2+Y^2+Z^2=1',
         'matrix_rule':'B0=T17^T*(G0+E00-I6/8)*T17+sum_event weight*row*row^T; B_alpha_sign=T^T*(I6/56+sign*eps^degree*E_alpha)*T/2'}
    (HERE/'six_root_general_formulas.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':
    export_general_formulas()
    verify_fixture()
