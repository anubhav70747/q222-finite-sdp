"""At most40 deterministic exact contact systems; no numerical optimization."""
from pathlib import Path
import json
import sympy as s

HERE=Path(__file__).resolve().parent
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
m,C,A,D,al,be=s.symbols('m C A D alpha beta')
unknowns=(m,C,A,D,al,be)
Du,Dv=1+u*u,1+v*v
F=(4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y*Y
 +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2)*Du*(v*Y-Z)**2
 +(4*be+2)*Dv*(X+u*Z)**2-2*Du*Dv*Z*Z)
eq=[s.diff(F,t)/2 for t in (W,X,Y,Z)]+[s.diff(F,u),s.diff(F,v)]
M,rhs=s.linear_eq_to_matrix(eq,unknowns)
evaluation=s.lambdify((u,v,W,X,Y,Z),(M,rhs),modules='sympy',cse=True)

def selector(cc,uc,vc,psi):
    mm,cc0,aa,dd,alpha,beta=cc
    dv=1+vc*vc; uu=uc*uc;du=1+uu
    ww,xx,yy,zz=psi
    av=s.Matrix([(ww+uc*yy)/du,(xx+uc*zz)/du])
    bv=s.Matrix([(ww-yy/uc)/du,(xx-zz/uc)/du])
    P=s.Matrix([[4*mm*vc*vc,-4*mm*vc],[-4*mm*vc,4*mm+(4*beta+2)*dv]])
    N=s.Matrix([[4*cc0*dv+4*dd+2,(4*dd+2)*vc],[(4*dd+2)*vc,(4*dd+2)*vc*vc]])
    rrhs=-P*av/uu; trhs=-uu*N*bv
    pp=s.cancel(rrhs[0]/(vc*vc*bv[0])); ss=s.cancel(trhs[1]/(vc*vc*av[1]))
    rr=s.cancel((rrhs[1]-ss*vc*vc*bv[1])/(vc*bv[0]))
    J=s.Matrix([[pp*vc*vc,0],[rr*vc,ss*vc*vc]])
    assert J*bv==rrhs and J.T*av==trhs
    aa0=4*mm;bb0=4*beta+2;cc1=4*cc0;ee=4*dd+2
    minor=s.factor(2*ee*cc1*(aa0+2*bb0)-ee*(rr-ss)**2-2*cc1*ss**2)
    return (pp,rr,ss),minor

results=[]
attempts=0
for hh in (s.Rational(1,2),s.Rational(2,3)):
    for uc in (s.Rational(1,2),s.Rational(1,3),s.Rational(1,4),s.Rational(1,5)):
        for vc in map(s.Integer,(2,3,4,6,8)):
            attempts+=1
            psi=(vc,-uc*hh,(1-hh)/vc,s.Integer(1))
            mat,rr=evaluation(uc,vc,*psi)
            if mat.det()==0:
                continue
            cc=tuple(mat.inv()*rr)
            if not all(z>0 for z in cc):
                continue
            pars,minor=selector(cc,uc,vc,psi)
            if minor>=0:
                continue
            score=sum(len(str(s.numer(q)))+len(str(s.denom(q))) for q in (*cc,*pars,minor))
            item={'u':str(uc),'v':str(vc),'state':list(map(str,psi)),
                  'source':dict(zip(map(str,unknowns),map(str,cc))),
                  'forced_selector':dict(zip(('p','r','s'),map(str,pars))),
                  'minor_at_v1':str(minor),'score':score}
            results.append(item)
assert attempts==40
results.sort(key=lambda a:a['score'])
output={'exact_systems_solved':attempts,'grid':'W=v,X=-u*h,Y=(1-h)/v; h in{1/2,2/3},u in{1/2,1/3,1/4,1/5},v in{2,3,4,6,8}',
        'negative_candidates':len(results),'candidates':results}
(HERE/'simple_source_candidates.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'exact_systems_solved':attempts,'negative_candidates':len(results),'best':results[:3]},indent=2))
