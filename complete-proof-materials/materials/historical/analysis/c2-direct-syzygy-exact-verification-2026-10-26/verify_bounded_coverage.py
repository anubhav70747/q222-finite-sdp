"""Two bounded additive checks from emitted data, no supplied-code imports."""
from pathlib import Path
import json,time
import sympy as s
START=time.perf_counter();OUT=Path(__file__).resolve().parent
d=json.loads((OUT/'private-copy/K364_SYZYGY_DATA.json').read_text())
u,v,W,X,Y,Z=s.symbols('u v W X Y Z',real=True)
n,m,C,A,D,alpha=s.symbols('n m C A D alpha',positive=True)
loc={str(x):x for x in (u,v,W,X,Y,Z,n,m,C,A,D,alpha)}
read=lambda x:s.sympify(x,locals=loc)
shape=s.Matrix([[read(x) for x in row] for row in d['base_row_shapes']]);k=[read(x) for x in d['base_row_factors']]
assert shape.shape==(16,4) and len(k)==16
full=s.diag(*k)*shape;full=full.col_join(v*full)
R=full[:,:3];q=full[:,3]
S=s.zeros(32,29)
for j,item in enumerate(d['syzygy_sparse']):
    for i,value in item['entries'].items():S[int(i),j]=read(value)
high=s.zeros(32,1)
for i,value in d['extra_top_degree_syzygy'].items():high[int(i)]=read(value)
lam=s.Matrix(s.symbols('ell0:29',real=True));tau=s.Symbol('tau',real=True)
t=s.diag(*[1/x for x in k+k])*(S*lam+tau*high)
assert all(s.simplify(s.expand(x))==0 for x in R.T*t)
# Expand the full completion with arbitrary independent t_i first. This
# leaves exactly the mixed-state syzygy term; no cancelling remainder is
# assumed after substituting the syzygy coordinates.
T=s.Matrix(s.symbols('T0:32',real=True));xi=s.Matrix([W,X,Y]);base=R*xi+q*Z
rho=-2*(1+u*u)*(1+v*v)**2
lhs=sum((base[i]+T[i]*Z)**2 for i in range(32))+(rho-2*(q.T*T)[0]-(T.T*T)[0])*Z**2
rhs=sum(base[i]**2 for i in range(32))+rho*Z**2
mixed=2*Z*(xi.T*R.T*T)[0]
assert s.expand(lhs-rhs-mixed)==0
# The separately verified R^T t=0 therefore kills precisely all three
# mixed-state coefficients for every source-only ell and tau.
print('PASS generic full four-state normal-completion identity: complete expansion leaves only 2 Z xi^T R^T t, and all three syzygy coefficients vanish',flush=True)
# Enumerate actual components of all 32 deformed roots in the fixed feature
# dictionary; no testing denominator, high-degree monomial, or omitted state.
features={(s.Rational(i,6),s.Rational(j,6),state) for i in range(7) for j in range(13) for state in range(4)}
root_terms=0;maxu=maxv=0;components=0
for i in range(32):
    for state in range(4):
        expr=R[i,state] if state<3 else q[i]+t[i]
        polynomial=s.Poly(s.expand(expr),u,v)
        components+=1
        for (eu,ev),coef in polynomial.terms():
            if coef==0:continue
            assert (s.Integer(eu),s.Integer(ev),state) in features,(i,state,eu,ev)
            assert not coef.has(u,v)
            root_terms+=1;maxu=max(maxu,eu);maxv=max(maxv,ev)
assert (maxu,maxv)==(1,2) and len(features)==364
# Pure-source contact-interpolation coefficients are a specialization of ell,
# and the zero-contact high correction still has the same represented form.
result={'result':'BOTH BOUNDED GENERIC COVERAGE CHECKS PASSED','elapsed_seconds':round(time.perf_counter()-START,6),'deformed_roots':32,'state_components_checked':components,'nonzero_component_monomials':root_terms,'maximum_testing_bidegree':[maxu,maxv],'fixed_library_features':len(features),'generic_source_only_gauge_coordinates':29,'extra_high_syzygy_parameter':True,'normal_completion_identity':'expanded exactly modulo the separately verified three-component syzygy','native_source_domain_or_contact_substitution_checked':False,'remaining_scalar_nonnegativity_or_SOS_checked':False}
(OUT/'bounded_coverage_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('PASS all 128 state components of 32 deformed roots use only fixed K364 dictionary features; maximum bidegree (1,2); arbitrary 29 source-only gauges and high parameter',flush=True)
print(json.dumps(result,indent=2),flush=True)
