"""Check the printed stationary dictionary against the original Bell slack.

Only small testing coefficients are cleared. The common source denominator
is removed before substitution; no feasibility solver or numerical test is used.
"""
from pathlib import Path
import contextlib, io, json, runpy
import sympy as s

HERE=Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    ctx=runpy.run_path(str(HERE/'full_r1_enlarged_gram.py'))
g=ctx['g']
for name in ['B','T','j','p','H','n','mu','m','C','A','D','alpha','fc','qc','rc','gc','kc','dc',
             'u','v','t','r','w','zeta','x','y','qz']:
    globals()[name]=ctx[name]
P=g['P']; source=g['src']; source_vars=list(source)
psource={a:s.cancel(P*b) for a,b in source.items()}
W=w-u*r; X=t-u*zeta; Y=r+u*w; Z=zeta+u*t; Du=1+u*u; Dv=1+v*v
half=s.Rational(1,2)
original=(m*(v*W-X+u*v*Y-u*Z)**2+n*Du*Dv*X*X+C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y
          +(D+half)*(u*W+u*v*X-Y-v*Z)**2+(alpha+half)*Du*(v*Y-Z)**2
          +half*Dv*(X+u*Z)**2-half*Du*Dv*Z*Z)
checks={}
def zero(name, expr):
    count=0
    for mon,co in s.Poly(s.expand(expr),u,v,t,r,w,zeta).terms():
        if co==0:
            continue
        reduced=0
        for powers,value in s.Poly(co,*source_vars).terms():
            assert sum(powers)<=1
            reduced+=value*(P if sum(powers)==0 else psource[source_vars[powers.index(1)]])
        assert s.cancel(reduced).as_numer_denom()[0]==0,(name,mon)
        count+=1
    checks[name]={'zero':True,'coefficients_cleared':count}

zero('original_Bell_slack_equals_complete_parent',original-Du*ctx['parent'])
Kd=C*H-(D+half)*(1-H); Kc=g['Kc']
for i,e in enumerate([
    T*(Kd+A*H+half)-A,
    j*Kd+mu+x*T*C*H,
    j*Kd+x*T*C*H+T*(C*H+A*H+(alpha+half)*H)-A-(alpha+half)*(1-p)+alpha/B,
    n*p-mu+p/2+(-mu-(alpha+half)*(1-p))*x+j*(alpha-n)+y*Kc]):
    zero('printed_scalar_cancellation_'+str(i+1),e)
(HERE/'source_parent_checks.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf8',newline='\n')
print('Original Bell slack and all four printed source cancellations pass exactly.')
