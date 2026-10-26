"""Only the graph coefficient definitions newly printed by refined fragment38."""
import json
from pathlib import Path
import sympy as s
u,v,W,Z,r,h,m,C,A,D,al,be=s.symbols('u v W Z r h m C A D al be')
X=-h*u*Z;Y=r*u*W;U=u*u;V=v*v;theta=1-h;half=s.Rational(1,2)
raw=m*(v*W-X+u*v*Y-u*Z)**2+C*(1+V)*(u*W-Y)**2+A*(1+U)*(1+V)*Y**2
raw+=(D+half)*(u*W+u*v*X-Y-v*Z)**2+(al+half)*(1+U)*(v*Y-Z)**2+(be+half)*(1+V)*(X+u*Z)**2-half*(1+U)*(1+V)*Z**2
Aw=((C+D+half)*(1-r)**2+A*r*r)*U+A*r*r*U*U
Aw+=(m+(2*m*r+C*(1-r)**2+(A+al+half)*r*r)*U+r*r*(m+A+al+half)*U*U)*V
Az=al+(al+(m+be+half)*theta**2)*U+(D+(2*D*h+half*h*h+be*theta**2)*U+(D+half)*h*h*U*U)*V
b0=m*theta+D*(1-r)+al*r+half
b1=m*r*theta+(D+half)*h*(1-r)+(al+half)*r
p=s.Poly(s.expand(raw-Aw*W**2-Az*Z**2+2*u*v*(b0+b1*U)*W*Z),W,Z)
out={'all_graph_coefficients_exact':all(cf==0 for cf in p.coeffs()),'scope':'Unconstrained original coefficients;all r,h;complete quarter-parent graph'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
