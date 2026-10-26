"""Bounded exact comparison with the recovered unequal-angle response text.

No source ZIP is used or claimed recovered. Requires SymPy and the frozen
independent general-angle JSON. Does not audit the order-six moment cone.
"""
import json
from pathlib import Path
import sympy as S

folder=Path(__file__).parent
prior=json.loads((folder/'c2_maxentangled_general_angles_2026_10_25.json').read_text())
a,r,s,z,c=S.symbols('a r s z c',positive=True)
e=r+s;p=r*s;g=a+e-1
L=a*(e-p)+p;Q=a*g+L
Delta=a**3*p+a*a*(e*e-p)+a*p*(2*e-p)+p*p
B=a**3*p+a*a*(e*e+p-e)+a*((3*e-2)*p-p*p)+p*p
source=[S.Integer(1),-a*B/(2*(a+1)*Delta),
        (a+s)*(1-r)*(a*s+r)*L/(2*(a+1)*g*Delta),
        r*r*(a+s)*Q/(2*(1-r)*Delta),
        a*(a+r)*(a+s)*L/(2*(a+1)*Delta),
        (a+r)*(1-s)*(a*r+s)*L/(2*(a+1)*g*Delta),
        s*s*(a+r)*Q/(2*(1-s)*Delta)]
local={'a':a,'r':r,'z':z}
old=[S.sympify(t,locals=local).subs(z,s/(a*r)) for t in prior['clean_source']]
assert all(S.cancel(u-v)==0 for u,v in zip(source,old))
helpers={name:S.sympify(t,locals=local).subs(z,s/(a*r))
         for name,t in prior['positive_helpers'].items()}
assert S.cancel(helpers['Q']-L/(a*r))==0
assert S.cancel(helpers['L']-g)==0
assert S.cancel(helpers['J']-(1-s))==0
assert S.cancel(helpers['R']-Q/a)==0
assert S.cancel(helpers['W']-Delta/(a*a*r*r))==0
assert S.cancel(helpers['N']+B/(a*a*r))==0
d,n,m,C,A,D,alpha=source
kappa=(a+1)*g/((1-r)*(1-s))
checks={
    'angle_ratio':a*(d+2*n)-2*A,
    'kappa_C':C+d/2-kappa*D,
    'kappa_alpha':alpha+d/2-kappa*m,
    'quadratic_relation':m*(C+d/2)-D*(alpha+d/2),
    'm_minus_D':m-D+(r-s)*L*Q/(2*(a+1)*g*Delta),
    'C_minus_alpha':C-alpha-(r-s)*L*Q/(2*(1-r)*(1-s)*Delta),
}
assert all(S.cancel(t)==0 for t in checks.values())
minor_old=S.sympify(prior['contact_minor_delete_d'],locals={'a':a,'c':c,'r':r})
Delta_nat=(a+1)*Delta/(a*a*r*r)
minor_new=-16384*c**6*r*r*(r-1)*(s-1)*g*g*Delta_nat/((1+c*c)**4*(1+a*a*c*c)**4)
assert S.cancel((minor_old-minor_new).subs(s,a*r*c*c))==0
witnesses=[]
for pair1,pair2 in [((S.Rational(1,8),S.Rational(1,16)),(S.Rational(1,16),S.Rational(1,8))),
                    ((S.Rational(1,16),S.Rational(1,8)),(S.Rational(1,8),S.Rational(1,16)))]:
    source1=[S.cancel(t.subs({a:1,r:pair1[0],s:pair1[1]})) for t in source]
    source2=[S.cancel(t.subs({a:2,r:pair2[0],s:pair2[1]})) for t in source]
    assert all(t>0 for t in source1) and all(t>0 for t in source2)
    R1=S.Matrix([[source1[2],source1[5]],
                  [source1[6]+source1[0]/2,source1[3]+source1[0]/2]])
    R2=S.Matrix([[source2[2],source2[5]],
                  [source2[6]+source2[0]/2,source2[3]+source2[0]/2]])
    k1=S.cancel(kappa.subs({a:1,r:pair1[0],s:pair1[1]}))
    k2=S.cancel(kappa.subs({a:2,r:pair2[0],s:pair2[1]}))
    determinant=(R1+R2).det()
    assert determinant==(k2-k1)*(source1[2]*source2[5]-source1[5]*source2[2])
    witnesses.append({'a1':1,'a2':2,'rs1':list(map(str,pair1)),'rs2':list(map(str,pair2)),
                      'source1':list(map(str,source1)),'source2':list(map(str,source2)),
                      'kappa1':str(k1),'kappa2':str(k2),
                      'mixture_determinant':str(determinant)})
expected=S.Rational(16557677,50331190)
assert sorted(S.Rational(t['mixture_determinant']) for t in witnesses)==[-expected,expected]
out={'input':'explicit recovered text only; no ZIP recovered',
     'source_order':['d','n','m','C','A','D','alpha'],
     'all_seven_symmetric_source_formulas_match':True,
     'six_helper_changes_verified':True,
     'source_domain_equivalence':'z=s/(a*r), old L=g, old N=-B0/(a^2*r), old J=1-s',
     'domain_includes_n_zero':'B0=0; no division by n or B0',
     'e_less_than_one_proof':'If e>=1 and 0<r,s<1 then all four terms in B0 are strictly positive: e^2+p-e>0 and (3e-2)-p>0. Contradiction.',
     'positive_denominators':'L>0 from e-p=r(1-s)+s>0; Q=ag+L>0; Delta0>0 termwise.',
     'angle_kappa_difference_identities_verified':{k:True for k in checks},
     'reciprocal_equivalence':'r=s iff m=D iff C=alpha; u*v=s/r, hence iff u*v=1',
     'extremality_minor_matches_independent_uncleared_operator_convention':True,
     'minor_nonzero_including_n_zero':'r,s<1,g>0,Delta_nat>0,c>0',
     'signed_witnesses':witnesses,
     'not_audited':'order-six moment pencil, degree-five matrix T, source ZIP, full conic hull'}
path=Path(__file__).with_suffix('.json')
path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k not in ['signed_witnesses']},indent=2))
print('Mixture determinants:',[t['mixture_determinant'] for t in witnesses])
