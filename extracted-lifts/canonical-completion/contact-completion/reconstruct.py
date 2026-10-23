"""Independent reconstruction from the completed Pro reply, without its archive.

Exact sparse rational arithmetic; no optimizer or numerical source tests.
Polynomial variable order: u,v,t,r,w,zeta,su,sv,sz.
Radical relations: su^2=u, sv^2=v, sz^2=w*zeta.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import json, hashlib

HERE=Path(__file__).resolve().parent
ZERO=(0,)*9
def add(*ps):
    out=defaultdict(F)
    for p in ps:
        for m,c in p.items():out[m]+=c
    return {m:c for m,c in out.items() if c}
def scale(p,c):return {m:x*c for m,x in p.items() if x*c}
def neg(p):return scale(p,-1)
def mul(p,q):
    out=defaultdict(F)
    for m,c in p.items():
        for n,d in q.items():
            e=[x+y for x,y in zip(m,n)]
            e[0]+=e[6]//2;e[6]%=2
            e[1]+=e[7]//2;e[7]%=2
            e[4]+=e[8]//2;e[5]+=e[8]//2;e[8]%=2
            out[tuple(e)]+=c*d
    return {m:c for m,c in out.items() if c}
def powp(p,n):
    out={ZERO:F(1)}
    for _ in range(n):out=mul(out,p)
    return out
def mon(*indices):
    e=[0]*9
    for i in indices:e[i]+=1
    return {tuple(e):F(1)}
one={ZERO:F(1)}
u,v,t,r,w,z,su,sv,sz=[mon(i) for i in range(9)]
def M(i,j):return mul(powp(u,i),powp(v,j))
W=add(w,neg(mul(u,r)));X=add(t,neg(mul(u,z)))
Y=add(r,mul(u,w));Z=add(z,mul(u,t))
Du=add(one,powp(u,2));Dv=add(one,powp(v,2))
H=[X,Y,Z,mul(v,W),mul(M(1,1),W),mul(u,W),mul(M(1,1),X),mul(v,Z),
   mul(M(2,1),W),mul(powp(u,3),W),mul(u,Z)]
strong=H+[mul(u,X),mul(v,X),mul(powp(v,2),W),mul(M(1,2),W),mul(v,Y),mul(M(1,1),Y),W,mul(M(1,1),Z),mul(u,Y)]
weak=[]
for root,exponents in [(t,[(0,0),(0,1),(1,0),(1,1)]),
                       (r,[(0,0),(0,1),(1,0),(1,1)]),
                       (w,[(0,0),(0,1),(1,0),(0,2),(1,1),(2,0),(0,3),(1,2),(2,1)]),
                       (z,[(0,0),(0,1),(1,0),(0,2),(1,1),(1,2),(2,1),(1,3),(2,2)])]:
    weak.extend(mul(M(i,j),root) for i,j in exponents)
weak.extend([mul(mul(su,sv),w),mul(mul(su,sv),z)])
weak.extend(mul(M(i,j),mul(su,sz)) for i,j in [(0,0),(0,1),(1,0),(1,1)])
weak.extend(mul(M(i,j),mul(sv,sz)) for i,j in [(0,0),(1,0)])
principal=[1,5,10,11,13,15,16,17,18,19,20,22,23,24,27,28]
strong_exp=[(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2),(3,0),(3,2),(4,0),(4,2),(5,0)]
weak_exp=[(1,0),(0,1),(2,0),(0,2),(1,1),(2,1),(2,2),(0,4)]
blocks=[('strong-base',strong,one)]
def gate(i,j):return scale(M(i,j),2**(12*i+25*j))
for i,j in strong_exp:
    for sign in [1,-1]:blocks.append((f'strong-{i}{j}-{sign:+}',H,add(one,scale(gate(i,j),sign))))
blocks.append(('weak-base',weak,Du))
for i,j in weak_exp:
    for sign in [1,-1]:blocks.append((f'weak-{i}{j}-{sign:+}',[weak[k-1] for k in principal],mul(Du,add(one,scale(gate(i,j),sign)))))

T={}
T['lambda']=mul(Du,powp(add(W,mul(v,X)),2))
T['n']=mul(mul(Du,Dv),powp(X,2))
T['m']=powp(add(mul(v,W),neg(X),mul(M(1,1),Y),neg(mul(u,Z))),2)
T['C']=mul(Dv,powp(add(mul(u,W),neg(Y)),2))
T['A']=mul(mul(Du,Dv),powp(Y,2))
T['D']=powp(add(mul(u,W),mul(M(1,1),X),neg(Y),neg(mul(v,Z))),2)
T['alpha']=mul(Du,powp(add(mul(v,Y),neg(Z)),2))
T['beta']=mul(Dv,powp(add(X,mul(u,Z)),2))
T['d']=scale(add(T['D'],T['alpha'],T['beta'],neg(mul(mul(Du,Dv),powp(Z,2)))),F(1,2))
sources=['d','lambda','n','m','C','A','D','alpha','beta']

columns=[]; allmon=set()
for b,(name,roots,weight) in enumerate(blocks,1):
    for i in range(len(roots)):
        for j in range(i,len(roots)):
            p=scale(mul(weight,mul(roots[i],roots[j])),1 if i==j else 2)
            columns.append((b,i+1,j+1,p));allmon.update(p)
for p in T.values():allmon.update(p)
monomials=sorted(allmon)
row_index={m:i+1 for i,m in enumerate(monomials)}
gram_records=[];source_records=[]
rows=defaultdict(list)
for b,i,j,p in columns:
    for monomial,c in sorted(p.items()):
        rec=[row_index[monomial],b,i,j,str(c)];gram_records.append(rec)
        rows[row_index[monomial]].append([b,i,j,str(c)])
for k in sources:
    for monomial,c in sorted(T[k].items()):source_records.append([row_index[monomial],k,str(c)])

# Construct literal integer-root positive-event Gram witnesses in strong block1.
event_roots={
 'lambda':[{18:1,13:1},{6:1,7:1}],
 'n':[{1:1},{12:1},{13:1},{7:1}],
 'm':[{4:1,1:-1,17:1,11:-1}],
 'C':[{6:1,2:-1},{5:1,16:-1}],
 'A':[{2:1},{20:1},{16:1},{17:1}],
 'D':[{6:1,7:1,2:-1,8:-1}],
 'alpha':[{16:1,3:-1},{17:1,11:-1}],
 'beta':[{1:1,11:1},{13:1,19:1}]
}
event_results={}
for name,rs in event_roots.items():
    polys=[add(*(scale(strong[k-1],c) for k,c in root.items())) for root in rs]
    event_results[name]=not add(*(powp(p,2) for p in polys),neg(T[name]))
assert all(event_results.values()),event_results

# Check all nine physical sign gauges in the ORIGINAL (unrotated) variables.
# Substitution is simultaneous, including the changed state signs.
def subst(p,images):
    out={}
    for exponents,coef in p.items():
        q={ZERO:coef}
        for k,n in enumerate(exponents):q=mul(q,powp(images[k],n))
        out=add(out,q)
    return out
# Build the original frame directly: w,t,r,z temporarily mean W,X,Y,Z.
OO={}
OO['lambda']=mul(Du,powp(add(w,mul(v,t)),2))
OO['n']=mul(mul(Du,Dv),powp(t,2))
OO['m']=powp(add(mul(v,w),neg(t),mul(M(1,1),r),neg(mul(u,z))),2)
OO['C']=mul(Dv,powp(add(mul(u,w),neg(r)),2))
OO['A']=mul(mul(Du,Dv),powp(r,2))
OO['D']=powp(add(mul(u,w),mul(M(1,1),t),neg(r),neg(mul(v,z))),2)
OO['alpha']=mul(Du,powp(add(mul(v,r),neg(z)),2))
OO['beta']=mul(Dv,powp(add(t,mul(u,z)),2))
OO['d']=scale(add(OO['D'],OO['alpha'],OO['beta'],neg(mul(mul(Du,Dv),powp(z,2)))),F(1,2))
gauge_results={}
for eps in [-1,1]:
    for eta in [-1,1]:
        images=[scale(u,eps),scale(v,eta),scale(t,eta),scale(r,eps),w,scale(z,eps*eta),su,sv,sz]
        gauge_results[f'{eps},{eta}']=all(not add(subst(p,images),neg(p)) for p in OO.values())

# Full-state sign-joining identity, cleared of all denominators.
# A=[[H,-uv],[-uv,H]]/2, B=[[-u²v,-uv²],[uv²,u²v]]/2.
# Verify 2 det(A)*K/Du after multiplying through explicitly via adjugate.
Hc=add(one,powp(u,2),powp(v,2));uv=M(1,1)
det4=add(powp(Hc,2),neg(powp(uv,2)))
hard=scale(add(mul(Hc,add(powp(t,2),powp(r,2))),scale(mul(uv,mul(t,r)),-2)),F(1,2))
cross=add(neg(mul(M(2,1),mul(t,w))),neg(mul(M(1,2),mul(t,z))),mul(M(1,2),mul(r,w)),mul(M(2,1),mul(r,z)))
B1=scale(add(neg(mul(M(2,1),w)),neg(mul(M(1,2),z))),F(1,2))
B2=scale(add(mul(M(1,2),w),mul(M(2,1),z)),F(1,2))
B_adj_B=scale(add(mul(Hc,add(powp(B1,2),powp(B2,2))),scale(mul(uv,mul(B1,B2)),2)),2)
scalar=scale(mul(mul(Du,Dv),add(mul(powp(uv,2),add(powp(w,2),powp(z,2))),scale(mul(mul(uv,Hc),mul(w,z)),-2))),F(1,2))
# K=Du*(hard+cross+B^T A^-1 B+scalar/det4).
sj5=add(mul(det4,T['d']),neg(mul(Du,mul(det4,add(hard,cross)))),neg(mul(Du,B_adj_B)),neg(mul(Du,scalar)))
assert not sj5

expected_samples={'t2':(t,25), 'r2':(r,17), 'w2':(w,9),'zeta2':(z,3)}
samples={}
for name,(root,expected_row) in expected_samples.items():
    exponent=next(iter(powp(root,2)))
    row=row_index[exponent]
    samples[name]={'row':row,'claimed_row':expected_row,'gram':rows[row],
                   'sources':[[k,str(T[k].get(exponent,0))] for k in sources if T[k].get(exponent,0)]}

def expected_diagonals(items):return sorted([[b,i,i,'1'] for b,i in items])
sample_expected={
 't2':(expected_diagonals([(b,1) for b in range(1,45)]),{'n':'1','m':'1','beta':'1','d':'1/2'}),
 'r2':(expected_diagonals([(b,2) for b in range(1,28)]+[(28,5)]+[(b,2) for b in range(29,45)]),{'C':'1','A':'1','D':'1','d':'1/2'}),
 'w2':(expected_diagonals([(1,18),(28,9)]),{'lambda':'1'}),
 'zeta2':(expected_diagonals([(b,3) for b in range(1,28)]+[(28,18)]+[(b,9) for b in range(29,45)]),{'alpha':'1'})
}
sample_results={name:(x['row']==x['claimed_row'] and sorted(x['gram'])==sample_expected[name][0]
                      and dict(x['sources'])==sample_expected[name][1]) for name,x in samples.items()}
assert all(sample_results.values()),sample_results
bounds=[1,1,512,8,16,4,16,128]
budget=sum(F(bound,2**(12*i+25*j)) for (i,j),bound in zip(weak_exp,bounds))
assert budget==F(2720379395100337099505665,9903520314283042199192993792)
B0=1-F(1,2**26);delta=F(1,2**26)
cap_containment=(1/B0**2-1 < F(1,2**24) and delta**2/(1-delta**2)<F(1,2**50))
assert cap_containment
# Strong baseline: use u temporarily as a separate source indeterminate h.
# Its only off-diagonal block is on rows4,6. After subtracting I/2,
# the determinant numerator with denominator4(1+h)^2 is17+14h-7h^2.
s=add(one,u)
strong_det_num=add(scale(mul(s,add(scale(one,4),neg(s))),7),scale(one,-4))
strong_det_identity=not add(strong_det_num,neg(add(scale(one,17),scale(u,14),scale(powp(u,2),-7))))
strong_fifth_identity=not add(add(scale(one,4),neg(powp(s,2))),
                              neg(mul(add(one,neg(u)),add(scale(one,3),u))))
assert strong_det_identity and strong_fifth_identity

result={'source_order':sources,'variable_order':['u','v','t','r','w','zeta','su','sv','sz'],
 'gram_block_orders':[len(r) for _,r,_ in blocks],
 'gram_order':sum(len(r) for _,r,_ in blocks),'total_order_with_signs':9+sum(len(r) for _,r,_ in blocks),
 'gram_variables':len(columns),'equations':len(monomials),
 'gram_records':len(gram_records),'source_records':len(source_records),
 'event_witnesses_exact':event_results,'all_nine_gauge_invariance':gauge_results,
 'full_state_sign_join_SJ5_exact':not sj5,'sample_rows':samples}
result.update({'all_four_sample_rows_exact':sample_results,'weak_budget_exact':str(budget),
               'weak_budget_below_1_over_1024':budget<F(1,1024),
               'remote_uncovered_cap_inside_contact_rectangle':cap_containment,
               'strong_baseline_minus_half_I_determinant_numerator_exact':strong_det_identity,
               'strong_baseline_fifth_diagonal_factor_exact':strong_fifth_identity,
               'strong_baseline_ge_half_I_for_0_le_H_le_1':True})
data={'conventions':'Rows, blocks and symmetric Gram indices are one-based. Source records are RHS.',
      'variable_order':result['variable_order'],'source_order':sources,
      'block_names':[n for n,_,_ in blocks], 'block_orders':result['gram_block_orders'],
      'monomials':[list(m) for m in monomials],
      'gram_records':gram_records,'source_records':source_records,'event_roots_block1':event_roots}
(HERE/'contact_affine_independent.json').write_text(json.dumps(data,separators=(',',':'))+'\n',encoding='utf8')
(HERE/'checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in result.items() if k!='sample_rows'},indent=2))
print('sample rows:',{n:(x['row'],x['claimed_row']) for n,x in samples.items()})
