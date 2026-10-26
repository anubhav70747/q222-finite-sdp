"""Exact SG7 contact check for two specified rational source families only."""
from pathlib import Path
import sympy as S
import json
s,x,w=S.symbols("s x w",positive=True)
q=S.Rational
def native(H,z,y,t,r):
    L=y*t; B=H+L; a=1+z*r; T=B+z*r; U=B-z*y*r; j=a-t; SS=a+L
    de=z+(1+z)*B+(1+z)*B*B+2*z*(1+z)*B*r+z**3*r*r+2*z*z*r+y*((1+z)*B+z+z*z*r*(2-r))
    gg=B+z*(H+r-1)
    mm=(1+z)*T*H+y*gg+y*z*j+z*r*L*(1+y)
    n=(T*a*SS-r*de)/(2*r*de)
    m=z*z*(1+y)*(1-r)*a*SS/(2*gg*de)
    C=a*mm/(2*L*de)
    A=y*a*T*SS/(2*H*de)
    D=(1+z)*T*U*SS/(2*y*j*de)
    al=T*mm/(2*(1-H)*de)
    return H,z,y,r,L,B,n,m,C,A,D,al
def check(vals):
    H,z,y,r,L,B,n,m,C,A,D,al=vals
    K=B+x*(H+r-1); J=L-w*(1+r*x)
    normals=[
      m*w*K+C*(1+w)*x*L+(D+q(1,2))*x*J,
      m*w*(1+x)*K+A*(1+x)*(1+w)*x*H+(al+q(1,2))*(1+x)*x*w*(H-1),
      -n*r*(1+x)*(1+w)*x*x*w-m*w*x*K+
            (D+q(1,2))*x*x*w*J+q(1,2)*(1+w)*x*x*w*(1-r)
    ]
    p=(n*(1+x)*(1+w)*r*r*x*x*w+m*w*K*K
      +C*(1+w)*x*L*L+A*(1+x)*(1+w)*H*H*x
      +(D+q(1,2))*x*J*J
      +(al+q(1,2))*(1+x)*x*w*(1-H)**2
      +q(1,2)*(1+w)*x*x*w*(1-r)**2
      -q(1,2)*(1+x)*(1+w)*x*w)
    exprs={"normal1_T0":normals[0],"normal2_T0":normals[1],
           "normal3_T0":normals[2],"parent_value":p,
           "partial_x":S.diff(p,x),"partial_w":S.diff(p,w)}
    result={name:str(S.cancel(expr.subs({x:z,w:y}))) for name,expr in exprs.items()}
    assert all(v=="0" for v in result.values()),result
    return result
out={
 "scope":"Exact stationary contacts plus nonnegative event weights; SG7 gives global validity.",
 "small_r_family":check(native(q(1,2),s,s,q(1,2),s)),
 "small_H_family":check(native(s,s*s,s*s,q(1,2),q(1,2))),
 "state_map_determinant":"plus or minus u^4*v^2; nonzero at positive source contact",
 "validity_reference":"latex/chapters/04_geometry.tex, SG7",
}
Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
print(json.dumps(out,indent=2))

