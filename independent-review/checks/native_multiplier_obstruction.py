"""Exact finite relabeling identities and the constants in the tight limit.
No maximization/solver: the existence of b_T is the compactness argument
in the matching note. Python and SymPy only.
"""
import json
from pathlib import Path
import sympy as S

def main():
    u,v,t=S.symbols('u v t')
    coeff=S.symbols('d lambda n m C A D alpha beta')
    def parent(cs):
        d,la,n,m,C,A,D,al,be=cs
        Du,Dv=1+u*u,1+v*v
        packets=[(4*la*Du,[1,v,0,0]),(4*n*Du*Dv,[0,1,0,0]),
                 (4*m,[v,-1,u*v,-u]),(4*C*Dv,[u,0,-1,0]),
                 (4*A*Du*Dv,[0,0,1,0]),(4*D+2*d,[u,u*v,-1,-v]),
                 ((4*al+2*d)*Du,[0,0,v,-1]),
                 ((4*be+2*d)*Dv,[0,1,0,u]),(-2*d*Du*Dv,[0,0,0,1])]
        out=S.zeros(4)
        for a,row in packets:
            z=S.Matrix(row);out+=a*z*z.T
        return out
    P=parent(coeff)
    events=[(0,1,1,1),(0,0,1,-1),(1,1,1,-1),(1,0,-1,1),
            (0,0,-1,1),(1,1,-1,1),(0,1,-1,-1),(1,0,1,-1)]
    maps=[lambda x,y,a,b:(1-x,y,a,(-1)**y*b),
          lambda x,y,a,b:(x,1-y,(-1)**(1-x)*a,b),
          lambda x,y,a,b:(y,x,(-1)**y*b,(-1)**x*a)]
    X=S.Matrix([[0,1],[1,0]])
    N=lambda w:S.Matrix([[1,w],[-w,1]])
    swap=S.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    data=[('gA',-u,-1/v,1+u*u,v*v,S.kronecker_product(N(u),S.eye(2))),
          ('gB',1/u,-v,1+v*v,u*u,S.kronecker_product(X,N(v))),
          ('gT',-1/v,-1/u,S.Integer(1),u*u*v*v,swap)]
    record={}
    cycle=[0,5,3,4,6,2,7,1]
    cycle_inverse={e:i for i,e in enumerate(cycle)}
    permutations=[]
    for func,(name,U,V,clear,factor,Nstate) in zip(maps,data):
        perm=[events.index(func(*e)) for e in events]
        permutations.append(tuple(perm))
        out=[coeff[0]]+[None]*8
        for i,j in enumerate(perm):out[j+1]=coeff[i+1]
        left=clear*parent(out)
        right=factor*Nstate.T*P.subs({u:U,v:V},simultaneous=True)*Nstate
        zero=(left-right).applyfunc(S.cancel)==S.zeros(4)
        assert zero
        cycleperm=[cycle_inverse[perm[e]] for e in cycle]
        record[name]={'raw_parent_identity':True,'event_permutation_zero_based':perm,
                      'event_cycle_permutation_zero_based':cycleperm}
    group={tuple(range(8))}
    while True:
        new=group|{tuple(g[h[i]] for i in range(8)) for g in group for h in permutations}
        if new==group:break
        group=new
    assert len(group)==16

    T=t*t
    source=[2*T-4*T**3,T**3,1-T+T**3,T**3,(1-T)/2+T**3,
            T**3,1-T+T**3,1-T+T**3,(1-T)/2+T**3]
    psi=S.Matrix([1+S.Rational(2,3)*t**4,-t**3/3,
                  -t+S.Rational(2,3)*t**3,-t**2-S.Rational(2,3)*t**4])
    test=S.expand((psi.T*parent(source).subs({u:-t,v:t})*psi)[0])
    poly=S.Poly(test,t)
    assert all(poly.coeff_monomial(t**i)==0 for i in range(6))
    assert poly.coeff_monomial(t**6)==S.Rational(8,3)
    norm=S.expand((1+t*t)**2*(psi.T*psi)[0])
    assert norm.subs(t,0)==1
    improvement=S.series(4*t**6-test/norm,t,0,8).removeO()
    assert S.expand(improvement).coeff(t,6)==S.Rational(4,3)
    M,kappa=S.symbols('M kappa')
    limit=S.Rational(4,3)*(M+2)-(4-kappa)
    assert S.expand(limit-S.Rational(4,3)*(M+3*kappa/4-1))==0
    result={'generator_identities':record,'finite_event_group_order':len(group),
            'cleared_state_pullback_t6_coefficient':str(poly.coeff_monomial(t**6)),
            'cleared_norm_constant':str(norm.subs(t,0)),
            'quantum_minus_local_t6_coefficient':str(improvement.coeff(t,6)),
            'tight_limit':str(S.factor(limit)),
            'kappa_limit_interval':['4/3','4'],
            'Motzkin_added_constant_interval':['0','2'],
            'scope':'Exact algebra supporting NATIVE_MULTIPLIER_OBSTRUCTION.md; the contact decomposition and limiting proof are analytical.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
