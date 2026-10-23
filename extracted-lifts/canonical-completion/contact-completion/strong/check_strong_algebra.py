"""Independent full-parent/division/conductor/shear/span audit.

The only pickle read is lII_generic.pkl freshly generated locally by the
inspected received extract_lII.py. No received pickle is loaded.
"""
from pathlib import Path
import pickle,json
import sympy as s
HERE=Path(__file__).resolve().parent
source=(HERE/'scripts/strong/extract_lII.py').read_text(encoding='utf8')
# Execute ONLY the already-inspected raw remainder definitions, before monic
# divisions or serialization. Do not rerun the generic extraction.
ctx={'__file__':str(HERE/'scripts/strong/extract_lII.py')}
exec(source.split('G=S.Poly(')[0],ctx)
for name in ['u','v','W','N','xi','ze','Rh','R','T','z','k','H','B','r','cp','dp','g','q','ss','ts','js','ls','als','al','nu','mu','aa','U','V']:
    globals()[name]=ctx[name]
d=pickle.loads((HERE/'data/lII_generic.pkl').read_bytes())
Q=d['Q'];Graw=ctx['G'];ff=s.Matrix([N,xi,ze,Rh,R,T])
Ractual=v*(u*u-z)*W;Tactual=u*k*(u*u-z)*W
Rhactual=(v*(u*xi-ze)+ss*u*(k*k-v*v)*W)/k
actual={R:Ractual,T:Tactual,Rh:Rhactual}
testing=(u,v,N,xi,ze,W)
results={}
def zero(expr,name,sub=None):
    terms=s.Poly(s.expand(expr),*testing)
    failures=[];count=0
    for exponents,c in terms.terms():
        if c==0:continue
        count+=1
        if sub:c=c.subs(sub,simultaneous=True)
        num,_=s.fraction(s.cancel(c))
        if num!=0:failures.append([list(exponents),str(s.factor(num))])
    results[name]={'coefficients_cleared':count,'zero':not failures,'failures':failures}
    print(name,count,'coefficients',len(failures),'failures',flush=True)
    assert not failures,(name,failures[:3])

# Exact generic missing-residual formula. dp,ts,js are absent from Q because
# their contributions occur only in these six contact remainder coefficients.
rx,rz,rn,ra,rb,rc=d['contact_residuals']
missing=(v*z*rx*xi*W+u*v*rz*ze*W+u*k*rn*N*W
         +(u**4*ra+u*u*rb+u*u*(v*v-k*k)*rc)*W*W)
qform=(ff.T*Q*ff)[0]
zero((Graw-qform-missing).subs(actual,simultaneous=True),'generic_full_reconstruction_with_all6_remainders')
assert not ({dp,ts,js}&Q.free_symbols)
results['Q_excludes_dp_ts_js_only_through_verified_contact_remainders']={'zero':True}

# Independently restore the complete original F and the four LII5 terms.
C0,D0,n0,m0,A0,alpha0=s.symbols('C0 D0 n0 m0 A0 alpha0')
eta=k*(N+q*Rh);X=xi-z*r*v*W/B;Y=eta+H*u*W/B;Z=ze+u*v*W/B
Du=1+u*u;Dv=1+v*v
J=u*W+u*v*X-Y-v*Z;M=v*W-X+u*v*Y-u*Z
full=(4*n0*Du*Dv*X*X+4*m0*M*M+4*C0*Dv*(u*W-Y)**2+4*A0*Du*Dv*Y*Y
      +4*D0*J*J+4*alpha0*Du*(v*Y-Z)**2
      +2*(J*J+Du*(v*Y-Z)**2+Dv*(X+u*Z)**2-Du*Dv*Z*Z))
retained=(4*(k*k*(C0+D0)-ls)*N*N+4*C0*v*v*eta*eta+4*n0*Du*Dv*xi*xi
          +4*(alpha0-als)*Du*(v*eta-ze)**2)
simple={cp:C0*k*k*ts,dp:D0*k*k*js,g:k*k*C0*D0/(C0+D0),q:D0/(C0+D0),
        ss:(1+z*r)/B,nu:n0*r,mu:m0/(z*z),aa:A0/(k*k),al:alpha0*(1-H)}
# Apply the pure-source linear relation B=H+y*t only after coefficient
# extraction; the testing/state polynomial remains small.
diff=s.expand((full-retained-Graw).subs(actual,simultaneous=True).subs(simple,simultaneous=True))
zero(diff,'original_F_equals_G_plus_all4_LII5_squares',{js:1+z*r-ts,B:H+k*k*ts})

# The added conductor is an exact pure-testing relation.
con=T*T-u*u*R*R-(u*u-z)/ss*Rh*T+u*(u*u-z)/ss*R*(u*xi-ze)
zero(con.subs(actual,simultaneous=True),'LII20_zero_conductor')

# Exact old/new source-only matrix congruence. u in its(2,3) entry is handled
# by the declared fixed testing span; all remaining shear factors are source.
hh=s.symbols('hh',positive=True)
new=s.Matrix([hh*N,xi+u*ze+R/(hh*hh*(1+hh*hh)),hh*hh*ze,hh*Rh,R/hh,T/hh])
mapback=s.diag(1/hh,1,1/(hh*hh),1/hh,hh,hh)
mapback[1,2]=-u/(hh*hh);mapback[1,4]=-1/(hh*(1+hh*hh))
assert all(s.cancel(x)==0 for x in mapback*new-ff)
Qc=s.hessian(con,list(ff))/2
assert s.expand((ff.T*Qc*ff)[0]-con)==0
results['exact_shear_map_and_conductor_hessian']={'zero':True}

# The actual fixed11/20 spans, in the original state variables.
Xs,Ys,Zs=s.symbols('Xs Ys Zs')
base=[Xs,Ys,Zs,v*W,u*v*W,u*W,u*v*Xs,v*Zs,u*u*v*W,u**3*W,u*Zs,
      u*Xs,v*Xs,v*v*W,u*v*v*W,v*Ys,u*v*Ys,W,u*v*Zs,u*Ys]
testspan=(u,v,W,Xs,Ys,Zs)
lookup={s.Poly(f,*testspan).monoms()[0]:i+1 for i,f in enumerate(base)}
eta0=Ys-H*u*W/B
xi0=Xs+z*r*v*W/B;ze0=Zs-u*v*W/B
Rh0=(u*v*Xs-v*Zs+ss*k*k*u*W)/k
N0=(Ys-q*u*v*Xs+q*v*Zs-(H/B+q*ss*k*k)*u*W)/k
newroots=[hh*N0,xi0+u*ze0+Ractual/(H*(1+H)),H*ze0,hh*Rh0,Ractual/hh,Tactual/hh]
span_results={}
def span(name,expr,principal=False):
    indices=[];bad=[]
    for monomial,c in s.Poly(s.expand(expr),*testspan).terms():
        if c==0:continue
        i=lookup.get(monomial)
        if i is None:bad.append(str(monomial))
        else:indices.append(i)
    ok=not bad and (not principal or max(indices)<=11)
    span_results[name]={'indices':sorted(set(indices)),'in_fixed_span':ok}
    assert ok,(name,span_results[name],bad)
for i,f in enumerate(newroots,1):span('normalized_root_'+str(i),f,True)
retained_roots={'N':N0,'v_eta':v*eta0,'xi':xi0,'u_xi':u*xi0,'v_xi':v*xi0,'uv_xi':u*v*xi0,
                'v_eta_minus_ze':v*eta0-ze0,'u_v_eta_minus_ze':u*(v*eta0-ze0)}
for name,f in retained_roots.items():span('retained_'+name,f)
results['root_spans']={'zero':True,'roots':span_results}
results['Q_free_symbols']=sorted(map(str,Q.free_symbols))
(HERE/'algebra_checks.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8')
print('All strong generic parent, conductor, shear and14 root-span checks passed.',flush=True)
