from pathlib import Path
import sympy as s
import json
p,C,A,al,b,z=s.symbols('p C A alpha b z',real=True)
W,Y,T=s.symbols('W Y T',real=True)
poly=4*C*(W-p*Y)**2+4*A*(1+p*p)*Y*Y+(4*al*(1+p*p)+2)*(Y-T)**2+2*p*p*Y*Y-4*p*W*T
M=s.hessian(poly,(W,Y,T))/2
subs={C:z/(2*(1-b*z)),A:(1+z*(1-b))/(2*b*z),al:1/(2*z*(1-b))}
result={'matrix':str(M),'minors':[str(s.factor(M[:i,:i].det())) for i in range(1,4)],
        'native_minors':[str(s.factor(M[:i,:i].det().subs(subs))) for i in range(1,4)]}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
