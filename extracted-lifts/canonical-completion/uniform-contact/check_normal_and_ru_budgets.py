"""Exact rational receipt for the independently derived R1 loss bounds."""
from fractions import Fraction as F
from pathlib import Path
import json

d=F(1,1024); S=F(81,65536); V=F(1,256); e=2*d; D=e*e+S*S
assert S>=d
f0=1-d-d**3; l0=1+d+d*d+d**3
R0=1+2*d+2*d*d+4*d**3+d**4*(1+2*d)
basic={
    'alpha_factor_lt_1':(1+d*d)*l0/(2*(1-d))<1,
    'c_factor_lt_1':l0/2<1,
    'fc_lower_1_8':f0/(4*R0)>F(1,8),
    'qc_lower_1_5':(1-2*d)/(4*R0)>F(1,5),
    'rc_lower_1_8':f0/(4*(1+d*d)*R0)>F(1,8),
    'phi_kap_fc_3T':F(8,3)*(1+d)*(1+3*d)/((1-d)*f0)<3,
    'phi_kap_sigma_24T':F(64,3)*(1+d)*(1+3*d)/(1-d)<24,
    'kap_qcy_3':(1+d)*(2+6*d)/((1-d)*(1-2*d))<3,
    'direct_phi_kap_rc_3T':F(4,3)*(2+6*d)/f0<3,
    'effective_Ru_lt_17':16+32*e*e<17,
}
rows={
    'Ru':[17*(32*R0/f0*S+F(512,3)*(d*d+d*S+S*S)),
          17*(64*d*S+F(256,3)*S*S),17*F(1280,3)*(e+S)**2,
          17*F(64,3)*R0*(1+d*d)/f0*S],
    'Jalpha':[16*(8*d*e*e+F(32,3)*e*e*S),F(0),F(0),16*F(64,3)*e*e*S],
    'O':[16*F(128,3)*d*V*(1+2*S)*D,16*32*V*V*D,
         16*F(512,3)*d*S*V*D,16*F(64,3)*d*V*D],
    'kc':[12*d+48*d*S,18*V,96*d*S,6*d],
    'H1_H2':[256*(d*d*V+V*V+d**3*S+d*S*V),F(0),256*(V*V+d*S*V),F(0)],
    'e0':[2*(d*d+V),F(0),F(0),F(0)],
}
totals=[sum(row[i] for row in rows.values()) for i in range(4)]
limits=[F(3,4),F(1,8),F(1,8),F(1,2)]
comparisons=[a<b for a,b in zip(totals,limits)]
# Polynomial identity proving the joint estimate, without sampling sources.
import sympy as s
p,q=s.symbols('p q'); j=p+q; J=j+1; K=J-p; Dp=1-p+p*p
diff=4*j*K*Dp-J*(j*(1-p)+2*p*p)
positive=(4*p*p-3*p+3)*q*q+(p+1)*(4*p*p-4*p+3)*q+3*p*(1-p)**2
assert s.expand(diff-positive)==0
assert s.expand(4*p*p-3*p+3-(4*(p-s.Rational(3,8))**2+s.Rational(39,16)))==0
assert s.expand(4*p*p-4*p+3-(4*(p-s.Rational(1,2))**2+2))==0
core=F(1,16)-d-F(9,4096)-V-S*V
assert core>F(1,32)
assert all(basic.values()) and all(comparisons)
result={'scope':'R1 source/purchase estimates only; no core, parent, support or coverage assertion',
        'basic_rational_checks':basic,'joint_source_polynomial_identity':True,
        'cap':{'u':'9/256','v':'1/16'},'core_reserve':str(core),
        'rows':{k:list(map(str,v)) for k,v in rows.items()},
        'totals':list(map(str,totals)),'limits':list(map(str,limits)),
        'strict_comparisons':comparisons}
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'basic_checks':basic,'strict_comparisons':comparisons,
                  'totals':list(map(str,totals)),'output':str(out)},indent=2))
