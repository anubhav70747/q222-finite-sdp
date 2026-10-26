"""Bounded exact reconstruction of sections 2--3 of the contact followup."""
import json
from pathlib import Path
import sympy as S

def main():
    a,b,c,r,hs,dd,nn,mm,CC,AA,DD,al=S.symbols('a b c r hs d n m C A D alpha', nonzero=True)
    delta=(b+c*c*r*hs)/a
    eps=S.Matrix([[0,-1],[1,0]])
    Psi=S.Matrix([[b/a,-c*r],[c*hs/a,1]])
    local=eps*Psi.T
    J=S.kronecker_product(S.eye(2),local.inv())
    basis=S.Matrix([[1,0,0,0],[0,1,0,-1],[0,1,0,1],[0,0,1,0]])
    invbasis=basis.inv()
    P22=S.diag(4*CC,4*nn+4*DD+2*dd,4*(mm+AA+al)+2*dd,0)
    # The leading raw coefficient is the physical matrix at A1=B1=-Z.
    Zi=S.diag(1,-1); Ii=S.eye(2)
    Aslim=[S.kronecker_product(Zi,Ii),S.kronecker_product(-Zi,Ii)]
    Bslim=[S.kronecker_product(Ii,Zi),S.kronecker_product(Ii,-Zi)]
    identity=S.eye(4)
    P22direct=dd*(identity-(Aslim[0]*Bslim[0]-Aslim[0]*Bslim[1]+Aslim[1]*Bslim[0]+Aslim[1]*Bslim[1])/2)
    evlim=[(1,-1,0,0),(1,-1,1,1),(-1,1,1,0),(-1,1,0,0),(-1,1,1,1),(-1,-1,0,1)]
    for weight,(sa,sb,i,j) in zip([nn,mm,CC,AA,DD,al],evlim):
        P22direct+=weight*(identity+sa*Aslim[i])*(identity+sb*Bslim[j])
    assert P22direct==P22
    K22=invbasis*J*P22*J.T*invbasis.T
    leading=S.factor(K22[2,3])
    claimed=-c*b*hs*(2*(mm+AA+al)+dd)/(a*a*delta*delta)
    assert S.factor(leading-claimed)==0

    ss,tt=S.symbols('s t')
    monomials=[1,ss,tt,ss*tt]
    unknown=S.symbols('v0:12')
    vv=[sum(unknown[4*i+j]*monomials[j] for j in range(4)) for i in range(3)]
    prod=S.Poly(vv[0]*ss+vv[1]*tt+vv[2]*ss*tt,ss,tt)
    coeff=S.Matrix([prod.coeff_monomial(ss**i*tt**j) for i,j in [(1,0),(0,1),(2,0),(0,2),(1,1),(2,1),(1,2),(2,2)]])
    linear=coeff.jacobian(unknown)
    syzygies=[[-tt,ss,0],[-tt,0,1],[-ss*tt,0,ss],[0,-ss*tt,tt]]
    syzcols=[]
    for row in syzygies:
        assert S.expand(row[0]*ss+row[1]*tt+row[2]*ss*tt)==0
        syzcols.append(S.Matrix([S.Poly(e,ss,tt).coeff_monomial(m) for e in row for m in monomials]))
    syzmat=S.Matrix.hstack(*syzcols)
    assert linear.rank()==8 and syzmat.rank()==4 and linear*syzmat==S.zeros(8,4)

    # Independently form the complete raw canonical parent for the rational witness.
    u,v=S.symbols('u v',real=True)
    I2=S.eye(2); Z=S.diag(1,-1); X=S.Matrix([[0,1],[1,0]])
    A=[Z,((1-u*u)*Z+2*u*X)/(1+u*u)]
    B=[Z,((1-v*v)*Z+2*v*X)/(1+v*v)]
    Ax=[S.kronecker_product(q,I2) for q in A]
    By=[S.kronecker_product(I2,q) for q in B]
    I4=S.eye(4)
    H=I4-(Ax[0]*By[0]-Ax[0]*By[1]+Ax[1]*By[0]+Ax[1]*By[1])/2
    weights=list(map(S.Rational,[0,'7/66','5/33','5/22','10/33','16/33','9/22',0]))
    events=[(1,1,0,1),(1,-1,0,0),(1,-1,1,1),(-1,1,1,0),(-1,1,0,0),(-1,1,1,1),(-1,-1,0,1),(1,-1,1,0)]
    physical=H
    for w,(sa,sb,i,j) in zip(weights,events):
        physical+=w*(I4+sa*Ax[i])*(I4+sb*By[j])
    P=physical.applyfunc(lambda q:S.cancel(q*(1+u*u)*(1+v*v)).expand())
    expected=P22.subs({dd:1,nn:weights[1],mm:weights[2],CC:weights[3],AA:weights[4],DD:weights[5],al:weights[6]})
    P22real=P.applyfunc(lambda q:S.Poly(q,u,v).coeff_monomial(u*u*v*v))
    assert P22real==expected
    Jw=J.subs({a:1,b:S.Rational(5,4),c:1,r:S.Rational(1,4),hs:S.Rational(1,2)})
    K=(invbasis*Jw*P*Jw.T*invbasis.T).applyfunc(S.expand)
    shifted=K.subs({u:ss+1,v:tt+1}).applyfunc(S.expand)
    zeta=shifted[:3,3]
    gamma=shifted[3,3]
    gp=S.Poly(gamma,ss,tt)
    G=S.Matrix([[gp.coeff_monomial(ss*ss),gp.coeff_monomial(ss*tt)/2,gp.coeff_monomial(ss*ss*tt)/2],
                [gp.coeff_monomial(ss*tt)/2,gp.coeff_monomial(tt*tt),gp.coeff_monomial(ss*tt*tt)/2],
                [gp.coeff_monomial(ss*ss*tt)/2,gp.coeff_monomial(ss*tt*tt)/2,gp.coeff_monomial(ss*ss*tt*tt)]])
    assert S.expand((S.Matrix([[ss,tt,ss*tt]])*G*S.Matrix([ss,tt,ss*tt]))[0]-gamma)==0
    leadmin=[G[:i,:i].det() for i in range(1,4)]
    assert all(q>0 for q in leadmin)
    def polar(mode):
        out=[]
        for pol in zeta:
            pp=S.Poly(pol,ss,tt)
            def co(i,j): return pp.coeff_monomial(ss**i*tt**j)
            assert co(0,0)==0
            if mode=='third_column':
                row=[co(1,0)+co(2,0)*ss,co(0,1)+co(0,2)*tt,co(1,1)+co(2,1)*ss+co(1,2)*tt+co(2,2)*ss*tt]
            elif mode=='equal_monomial_division':
                row=[co(1,0)+co(2,0)*ss+co(1,1)*tt/3+co(2,1)*ss*tt/2,
                     co(0,1)+co(0,2)*tt+co(1,1)*ss/3+co(1,2)*ss*tt/2,
                     co(1,1)/3+co(2,1)*ss/2+co(1,2)*tt/2+co(2,2)*ss*tt]
            elif mode=='first_two_st':
                row=[co(1,0)+co(2,0)*ss+co(1,1)*tt/2+co(2,1)*ss*tt/2,
                     co(0,1)+co(0,2)*tt+co(1,1)*ss/2+co(1,2)*ss*tt/2,
                     co(2,1)*ss/2+co(1,2)*tt/2+co(2,2)*ss*tt]
            out.append(row)
        out=S.Matrix(out)
        assert (out*S.Matrix([ss,tt,ss*tt])-zeta).applyfunc(S.expand)==S.zeros(3,1)
        return out
    witness={}
    for mode in ['third_column','equal_monomial_division','first_two_st']:
        V=polar(mode).subs({ss:-1,tt:0})
        remainder=K[:3,:3].subs({u:0,v:1})-V*G.inv()*V.T
        witness[mode]=str(S.factor(remainder[0,0]))
    result={'leading_zeta3':str(leading),'leading_zeta3_claim_verified':True,'P22_verified_generically':True,'P22_verified_on_witness':True,
            'bilinear_coefficient_map_rank':8,'bilinear_syzygy_dimension':4,
            'claimed_syzygies_basis':True,'G_witness':[[str(q) for q in G.row(i)] for i in range(3)],
            'G_witness_leading_principal_minors':list(map(str,leadmin)),
            'canonical_polarization_candidate_values':witness,
            'claimed_negative_value':'-6975688592/51747662835'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
