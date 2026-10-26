"""Small exact algebra for the half-power contact normal form.
No parameter test, positivity solver, or universal admission claim.
"""
import json
from pathlib import Path
import sympy as S

def main():
    s,t,uc,vc=S.symbols('s t uc vc')
    q=S.Matrix([s,t,s*t]);u=uc+s;v=vc+t
    Mu=S.Matrix([[0,0,0],[0,uc,S.Rational(1,2)],[0,S.Rational(1,2),0]])
    Mv=S.Matrix([[vc,0,S.Rational(1,2)],[0,0,0],[S.Rational(1,2),0,0]])
    assert S.expand((q.T*Mu*q)[0]-u*t*t)==0
    assert S.expand((q.T*Mv*q)[0]-v*s*s)==0
    coeff=S.symbols('c0:16');mon=[1,s,t,s*t]
    V=[sum(coeff[4*i+j]*mon[j] for j in range(4)) for i in range(3)]
    ell=coeff[12]+coeff[13]*t;mm=coeff[14]+coeff[15]*s
    pol=S.Poly(sum(V[i]*q[i] for i in range(3))+u*t*ell+v*s*mm,s,t)
    exps=[(1,0),(0,1),(2,0),(0,2),(1,1),(2,1),(1,2),(2,2)]
    matrix=S.Matrix([pol.coeff_monomial(s**i*t**j) for i,j in exps]).jacobian(coeff)
    assert matrix.rank()==8
    assert S.expand((-u*ell)*t+u*t*ell)==0
    assert S.expand((-v*mm)*s+v*s*mm)==0
    parity={(a,b):[] for a in range(2) for b in range(2)}
    full={(a,b):[] for a in range(2) for b in range(2)}
    for state in 'WXYZ':
        for i in range(3):
            for j in range(3):
                key=(i%2,j%2)
                full[key].append((state,i,j))
                if state=='W' and i+j<2:continue
                if state=='Z' and i+j>2:continue
                parity[key].append((state,i,j))
    order=[(0,0),(1,0),(0,1),(1,1)]
    raw=[len(full[x]) for x in order]
    pruned=[len(parity[x]) for x in order]
    contact=[n-1 for n in pruned]
    assert raw==[16,8,8,4] and pruned==[14,6,6,4] and contact==[13,5,5,3]
    result={'scalar_Mu_identity':True,'scalar_Mv_identity':True,
            'cross_coefficient_unknowns_per_row':16,'cross_map_rank':8,
            'free_cross_parameters_per_row':8,'total_free_cross_parameters':24,
            'new_weighted_channel_cancellations':True,
            'fixed_unshifted_parity_orders':raw,'fixed_total':sum(raw),
            'contact_before_support_orders':[15,7,7,3],
            'C2_support_pruned_orders':pruned,
            'C2_support_and_contact_orders':contact,
            'C2_support_and_contact_total':sum(contact),
            'literal_kept_root_monomials':{str(k):parity[k] for k in order},
            'scope':'Exact normal-form identities/counts only. Universal PSD Gram feasibility is not tested or asserted.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='literal_kept_root_monomials'},indent=2))
if __name__=='__main__':main()
