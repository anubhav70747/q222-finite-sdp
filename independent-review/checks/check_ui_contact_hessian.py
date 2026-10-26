"""Independent matched UI Hessian from RC4. No numerical sampling."""
import json
import sympy as s

a,b,r=s.symbols('a b r',nonzero=True)
u,V,X,Z=s.symbols('u V X Z')
n,m,al,H,Aw=s.symbols('n m al H Aw')
z=a*a; Ac=1+z*r; g=(1+z)*b-1+r
De=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r
sources={n:((b+r)*Ac*Ac-r*De)/(2*r*De),m:(1-r)*Ac*Ac/(2*g*De),
         al:z*z*(1+z)*b*(b+r)**2/(2*(1-z*b)*De),
         H:z*(1+z)*(b+r)*b/(2*De),Aw:Ac*Ac*(b+r)/(2*z*b*De)}
f=(4*n*(1+u*u)*X*X+4*m*(V*(1+u*u)-X-u*Z)**2
   +4*H*V*V*(u*X-Z)**2+4*Aw*(1+u*u)*u*u
   +(4*al+2)*(1+u*u)*(V*u-Z)**2
   +2*(X+u*Z)**2-2*(1+u*u)*Z*Z)
q0={u:a,V:1,X:-r/b,Z:1/(a*b)}
hh=(s.hessian(f,(X,Z,V,u)).subs(q0)/2).applyfunc(lambda v:s.cancel(v.subs(sources)))
print('source Hessian reconstructed',flush=True)
field=s.QQ.frac_field(a,b,r)
ff=field.from_sympy
HH=[[ff(hh[i,j]) for j in range(4)] for i in range(4)]
dp=HH[0][0]*HH[1][1]-HH[0][1]*HH[1][0]
adj=[[HH[1][1],-HH[0][1]],[-HH[1][0],HH[0][0]]]
expected=s.Matrix([
 [8*Ac*(z*z*b*b+z*z*b*r+z*b*b+2*z*b*r-z*b-r+1)/(b*De),
  4*Ac*Ac*(z*b+z*r+2*r-2)/(a*b*De)],
 [4*Ac*Ac*(z*b+z*r+2*r-2)/(a*b*De),
  8*(1-r)*Ac*Ac/(z*b*De)]])
checks={}
for i,j in ((0,0),(0,1),(1,1)):
    correction=sum((HH[l][i+2]*adj[l][ll]*HH[ll][j+2] for l in range(2) for ll in range(2)),field.zero)
    checks[f'UI_G_{i}{j}']=(HH[i+2][j+2]-ff(expected[i,j]))*dp-correction==field.zero
    print(f'checked {i}{j}',flush=True)
det_expected=16*(b+r)**2*Ac**3*(3*z+4-(z+2)**2*r)/(b*b*De*De)
checks['UI5_determinant']=ff(expected[0,0])*ff(expected[1,1])-ff(expected[0,1])**2-ff(det_expected)==field.zero
Qn=1+r*(z-1)-z*(z+2)*r*r-z*(z+1)*b*r
checks['Qn_numerator']=s.cancel((b+r)*Ac*Ac-r*De-b*Qn)==0
checks['Qn_factorization']=s.expand(Qn-((1-(z+1)*r)*Ac-z*r*g))==0
checks['UI6_endpoint_bound']=s.cancel(3*z+4-(z+2)**2/(z+1)-z*(2*z+3)/(z+1))==0
checks['gradient_contact']=[s.cancel(s.diff(f,t).subs(q0).subs(sources))==0 for t in (X,Z,V,u)]
checks['value_contact']=s.cancel(f.subs(q0).subs(sources))==0
out={'checks':checks,'all_pass':all(all(v) if isinstance(v,list) else v for v in checks.values()),
     'scope':'RC4 literal leading polynomial with each testing a replaced by u; no scalar/normal allocation replay.'}
print(json.dumps(out,indent=2))

raise SystemExit(0 if out["all_pass"] else 1)
