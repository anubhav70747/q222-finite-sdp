"""Exact zero that excludes the earlier separated Alice-arc construction."""
from pathlib import Path
import json,runpy
import sympy as s
HERE=Path(__file__).resolve().parent
g=runpy.run_path(str(HERE/'source_dictionary.py'))
f=s.Rational
par={g['B']:f(1,2**20),g['T']:f(1,2**40),g['j']:1,g['p']:f(1,2),g['H']:f(1,2)}
source={str(a):s.factor(b.subs(par)) for a,b in g['src'].items()}
assert all(source[k]>0 for k in ['n','m','C','A','D','alpha'])
assert s.expand(g['NN'].subs({g['j']:1,g['p']:f(1,2),g['H']:f(1,2)})
                -(f(1,2)-g['B']/4-g['B']**2*g['T']/4))==0
u=f(1,2**10);v=f(1,2**30);Du=1+u*u;Dv=1+v*v
W=1+f(1,2**61);X=-f(1,2**31)
Y=f(1,2**10)*(1-f(1,2**41));Z=f(1,2**20)*(1+f(1,2**21))
original=(source['m']*(v*W-X+u*v*Y-u*Z)**2+source['n']*Du*Dv*X*X
          +source['C']*Dv*(u*W-Y)**2+source['A']*Du*Dv*Y*Y
          +(source['D']+f(1,2))*(u*W+u*v*X-Y-v*Z)**2
          +(source['alpha']+f(1,2))*Du*(v*Y-Z)**2
          +f(1,2)*Dv*(X+u*Z)**2-f(1,2)*Du*Dv*Z*Z)
assert original==0
assert f(1,2**24)<u*u<f(2**24) and 0<v<f(1,2**25)
assert abs(u)<f(9,256) and abs(v)<f(1,16)
nX2=source['n']*X*X
C_cE2=source['C']*(v*v/Dv)*(u*W-Y)**2/Du
assert nX2>0 and C_cE2>0
assert (1-f(1,2**26))**2*Du>1 # b=1/sqrt(Du)<1-2^-26
assert v*v/Dv<f(1,2**52) # Bob-complement threshold is not reached
out={'source_parameters':{str(a):str(b) for a,b in par.items()},
     'state':list(map(str,[W,X,Y,Z])),'u':str(u),'v':str(v),
     'source_weights':{a:str(b) for a,b in source.items()},
     'original_parent_is_zero':True,'inside_earlier_testing_annulus':True,
     'inside_enlarged_contact_cap':True,
     'retained_nX_squared':str(nX2),'retained_CcE_squared':str(C_cE2),
     'scope':'The earlier separated Alice-arc Gram cannot admit this source; other completions are not excluded.'}
(HERE/'annulus_obstruction.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
print('Exact annulus zero and both strictly positive retained terms verified.')
