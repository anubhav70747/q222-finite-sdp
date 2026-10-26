"""One-source exact facial reduction and rational witness recovery.

No source search.  The source is the fixed scaled simple C1 example.
"""
from pathlib import Path
import json
import hashlib
import sympy as s

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json').read_text())
diag=json.loads((HERE/'simple_source_numerical_diagnostic.json').read_text())
names=data['auxiliary_order']
y=s.symbols(' '.join(names))
source=list(map(s.Rational, ['359','0','0','75/16','60','416','41/2','831/2','151/2']))
params=source+list(y)
full=s.zeros(24)
for k,i,j,a in data['entries']:
    full[i,j]+=a*params[k]
SG=full[:8,:8];SH=full[8:16,8:16];TG=full[16:20,16:20];TH=full[20:24,20:24]
g=s.Matrix([34,4,6,-9]);h=s.Matrix([17,2,12,-18])
eq=list(SG[0,:])+list(SG[3,:])
eq+=list(SG*(g.col_join(2*g)))+list(TG*g)
eq+=list(SH*(h.col_join(2*h)))+list(TH*h)
# SG's endpoint zero rows imply TG has zero diagonal entries 0 and 3.
# Positivity then forces these two entire rows to vanish as well.
eq+=list(TG[0,:])+list(TG[3,:])
A,b=s.linear_eq_to_matrix(eq,y)
solset=s.linsolve((A,b),y)
if solset is s.EmptySet: raise AssertionError('inconsistent facial equations')
sol=list(next(iter(solset)))
free=sorted(set().union(*(a.free_symbols for a in sol)),key=lambda a:names.index(str(a)))
print('Face equations:',A.rows,'rank:',A.rank(),'free:',len(free))
print('Free coordinates:',list(map(str,free)))
subs=dict(zip(y,sol))
blocks=[M.subs(subs,simultaneous=True) for M in (SG,SH,TG,TH)]
target={str(a):359*v for a,v in zip(y,diag['auxiliary'])}
trial={a:s.Rational(str(target[str(a)])).limit_denominator(64) for a in free}
print('One nearest-rational trial:',{str(a):str(b) for a,b in trial.items()})
def psd_ldl(M):
    """Exact symmetric elimination; returns failure or a rational LDL congruence."""
    n=M.rows
    rest=s.Matrix(M);L=s.eye(n);ds=[]
    for i in range(n):
        p=s.factor(rest[i,i]);ds.append(p)
        if p<0:return {'ok':False,'index':i,'pivot':str(p)}
        if p==0:
            if any(rest[i,j]!=0 for j in range(i+1,n)):
                return {'ok':False,'index':i,'zero_pivot_nonzero_row':True}
            continue
        for j in range(i+1,n):L[j,i]=s.cancel(rest[j,i]/p)
        for j in range(i+1,n):
            for k in range(j,n):
                rest[j,k]=s.cancel(rest[j,k]-rest[j,i]*rest[i,k]/p)
                rest[k,j]=rest[j,k]
    assert M==L*s.diag(*ds)*L.T
    return {'ok':True,'L':L,'D':ds}
candidate=[M.subs(trial) for M in blocks]
certs=[psd_ldl(M) for M in candidate]
print('One rational trial PSD:',[{k:v for k,v in z.items() if k not in ['L','D']} for z in certs])
literal=[s.cancel(a.subs(trial)) for a in sol]
print('Selector:',list(map(str,literal[:12])))
print('Ranks:',[sum(bool(a) for a in z['D']) if z['ok'] else None for z in certs])
print('Positive pivots:',[[float(a) for a in z['D'] if a>0] if z['ok'] else None for z in certs])
# A second exact coordinate description puts selector coordinates last so
# rounding does not inflate their displayed denominators.  This is one
# deterministic rounding operation on the already recovered witness, not
# a new source or a numerical conic solve.
ordered=list(y[12:])+list(y[:12])
AA,bb=s.linear_eq_to_matrix(eq,ordered)
ss=list(next(iter(s.linsolve((AA,bb),ordered))))
sf=sorted(set().union(*(a.free_symbols for a in ss)),key=lambda a:ordered.index(a))
old=dict(zip(y,literal))
grid=lambda value,N:s.Rational(s.floor(N*value+s.Rational(1,2)),N)
st={a:grid(old[a],64) for a in sf}
newdict={a:s.cancel(b.subs(st)) for a,b in zip(ordered,ss)}
newliteral=[newdict[a] for a in y]
newblocks=[M.subs(newdict,simultaneous=True) for M in (SG,SH,TG,TH)]
newcerts=[psd_ldl(M) for M in newblocks]
assert all(z['ok'] for z in newcerts)
print('Uniform-grid selector:',list(map(str,newliteral[:12])))
print('Uniform-grid free:',{str(a):str(b) for a,b in st.items()})
print('Uniform-grid ranks:',[sum(bool(a) for a in z['D']) for z in newcerts])
payload={'source':list(map(str,source)),'source_scale_relative_d1':'359',
    'auxiliary_order':names,'face_equation_count':A.rows,'face_rank':A.rank(),
    'free_coordinates':list(map(str,free)),
    'affine_solution':list(map(str,sol)),
    'first_rational_trial':{str(a):str(b) for a,b in trial.items()},
    'first_trial_psd':[bool(z['ok']) for z in certs]}
(HERE/'simple_source_face_reduction.json').write_text(json.dumps(payload,indent=2)+'\n')

def serial(M):return [[str(a) for a in row] for row in M.tolist()]
def compact(M):
    den=s.ilcm(*(s.denom(a) for a in M))
    return {'denominator':str(den),'integer_numerator':serial(den*M)}
witness={'format':'q222-fixed-simple-source-halfline-rational-witness-v1',
 'scope':'One fixed valid C1 source. No claim that quadratic half-line selectors suffice for every valid source.',
 'source_order':data['source_order'],'source':list(map(str,source)),
 'source_d1':list(map(str,[a/359 for a in source])),
 'auxiliary_order':names,'auxiliary':list(map(str,newliteral)),
 'selector_common_denominator':'64',
 'selector_integer_coefficients':[[[0,0],[1403,0]],[[3043,0],[23874,815]],[[-2109,3875],[2306,43455]]],
 'BG':serial(s.Matrix(4,4,newliteral[12:28])),
 'BH':serial(s.Matrix(4,4,newliteral[28:44])),
 'BG_compact':compact(s.Matrix(4,4,newliteral[12:28])),
 'BH_compact':compact(s.Matrix(4,4,newliteral[28:44])),
 'blocks':{},'ldl':{},'face_rank':30,'face_dimension':14,
 'original_pencil_sha256':hashlib.sha256((HERE/'HALFLINE_QUADRATIC_SELECTOR_ORDER24.coo.json').read_bytes()).hexdigest()}
for name,M,z in zip(('SG','SH','TG','TH'),newblocks,newcerts):
    witness['blocks'][name]=serial(M)
    witness['ldl'][name]={'L':serial(z['L']),'D':list(map(str,z['D'])),
                         'rank':sum(bool(a) for a in z['D']),
                         'positive_pivot_floor':[int(s.floor(a)) for a in z['D'] if a>0]}
(HERE/'SIMPLE_SOURCE_HALFLINE_EXACT_WITNESS.json').write_text(json.dumps(witness,indent=2)+'\n')
print('Compact BG:',witness['BG_compact'])
print('Compact BH:',witness['BH_compact'])
print('Positive pivot floors:',{a:z['positive_pivot_floor'] for a,z in witness['ldl'].items()})
