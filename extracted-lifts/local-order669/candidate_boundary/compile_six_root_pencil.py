"""Evaluate the exact affine Gram pencil for the fixed six-root testing family.

No SDP is solved. This cone certifies nonnegativity on the stated LOCAL test
rectangle. It is not asserted to be a sound cone of globally valid Bell slacks.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import sys
from collections import defaultdict
import compile_six_root_chart as chart

S=chart.S
HERE=Path(__file__).resolve().parent
u,w,W,X,Y,Z=S.symbols('u w W X Y Z')
states=(W,X,Y,Z)
generators=(u,w,W,X,Y,Z)
env=dict(zip(('u','w','W','X','Y','Z'),generators))
features=[tuple(S.Poly(S.sympify(expr,locals=env),*generators).monoms()[0]) for expr in chart.PHI]
structured='--structured' in sys.argv
prefix='six_root_structured' if structured else 'six_root'
feature_blocks=([features[:10]]*15+[[features[i] for i in inds] for inds,_ in chart.EVENT_GROUPS]) if structured else ([features]+[features[:10]]*14)
weight_blocks=[{(0,0):F(1)}]
for alpha in chart.ALPHAS[1:]:
    for sign in (1,-1):
        weight_blocks.append({(0,0):F(1),alpha:F(sign*4096**sum(alpha))})
if structured:weight_blocks += [{(0,0):F(1)}]*4
sizes=[len(v) for v in feature_blocks]
offsets=[];offset=0
for n in sizes:offsets.append(offset);offset+=n

positions=[]
rows=defaultdict(lambda:defaultdict(F))
for block,(flist,weights) in enumerate(zip(feature_blocks,weight_blocks)):
    for i,fi in enumerate(flist):
        for j in range(i,len(flist)):
            fj=flist[j]; col=len(positions);positions.append((block,i,j))
            p=tuple(a+b for a,b in zip(fi,fj))
            mult=1 if i==j else 2
            for alpha,coef in weights.items():
                e=(p[0]+alpha[0],p[1]+alpha[1])+p[2:]
                rows[e][col]+=coef*mult
count=len(positions)
cs=S.symbols('d m C A D alpha beta')
d,m,C,A,D,ap,beta=cs
Du=1+u*u;Dw=1+w*w
parent=4*m*(W-w*X+u*Y-u*w*Z)**2+4*C*Dw*(u*W-Y)**2+4*A*Du*Dw*Y*Y+(4*D+2*d)*(u*w*W+u*X-w*Y-Z)**2+(4*ap+2*d)*Du*(Y-w*Z)**2+(4*beta+2*d)*Dw*(X+u*Z)**2-2*d*Du*Dw*Z*Z
target={}
for k,c in enumerate(cs):
    for mon,coef in S.Poly(S.diff(parent,c),*generators).terms():
        val=F(str(coef));rows[mon][count+k]-=val
        target.setdefault(mon,{})[k]=val

def add_scaled(dest,source,mult):
    for k,v in source.items():
        val=dest.get(k,F(0))+mult*v
        if val:dest[k]=val
        elif k in dest:del dest[k]

original={mon:dict(vals) for mon,vals in rows.items()}
pivots={}
for mon in sorted(rows):
    row={k:v for k,v in rows[mon].items() if v}
    while row:
        p=min(row)
        if p>=count:
            raise AssertionError(('target has unresolved visible constraint',mon,row))
        if p in pivots:
            add_scaled(row,pivots[p],-row[p])
        else:
            scale=row[p]
            pivots[p]={k:v/scale for k,v in row.items()}
            break
free=[k for k in range(count) if k not in pivots]
free_index={k:j for j,k in enumerate(free)}
# Affine coefficients use1..7 for source coordinates,8+j for auxiliaries.
expressions={k:{8+j:F(1)} for k,j in free_index.items()}
for pivot in sorted(pivots,reverse=True):
    ans={}
    for k,v in pivots[pivot].items():
        if k==pivot:continue
        if k>=count:add_scaled(ans,{1+k-count:F(1)},-v)
        else:add_scaled(ans,expressions[k],-v)
    expressions[pivot]=ans

# Independently reconstruct every polynomial coefficient from the final pencil.
for mon,row in original.items():
    residue={}
    for k,v in row.items():
        if k<count:add_scaled(residue,expressions[k],v)
        else:add_scaled(residue,{1+k-count:F(1)},v)
    assert not residue,(mon,residue)

entries=[]
for col,(block,i,j) in enumerate(positions):
    for coeff,val in sorted(expressions[col].items()):
        if val:entries.append([coeff,offsets[block]+i,offsets[block]+j,str(val)])
output={'status':'Exact LOCAL Gram cone; nonnegativity only on0<=u,w<=1/4096; no global Bell-cone or Q222 equality asserted',
        'canonical_source_order':[str(c) for c in cs],
        'matrix_order':offset,'auxiliary_count':len(free),'matching_rank':len(pivots),
        'raw_gram_variable_count':count,'coefficient_equation_count':len(original),
        'block_orders':sizes,'blocks':[{'index':i,'offset':o,'order':n} for i,(o,n) in enumerate(zip(offsets,sizes))],
        'entries':entries,
        'coefficient_convention':'[k,i,j,rational], upper triangle, mirror once; k0constant,k1..7canonical source,k8+jauxiliaryj',
        'feature_family':'six_root_general_formulas.json; base17 then +/- for each listed nonzero monomial',
        'normalization':'unnormalized reciprocal parent; multiply allfeatures by kappa_t to obtain normalized Bell slack',
        'positive_domain':'0<=u,w<=1/4096, all realWXYZ; normalized state sphere optional',
        'free_gram_positions':[list(positions[i]) for i in free],
        'source_admission_proved':'nativeP1guards SR4-SR5; inparticular SR7 x,sigma<=1/1024,z<=1,lambda<=0'}
if structured:
    output['feature_family']='15principal10functions;fourretainedgroups incompile_six_root_chart.EVENT_GROUPS'
    output['per_block_feature_monomials']=[[list(m) for m in f] for f in feature_blocks]
(HERE/(prefix+'_source_pencil.coo.json')).write_text(json.dumps(output,indent=2)+'\n')
product={'monomial_order':'u,w,W,X,Y,Z','gram_positions':[list(x) for x in positions],
         'coefficient_rows':[{'monomial':list(mon),'gram':[[k,str(v)] for k,v in sorted(row.items()) if k<count],
                             'target':[[k-count,str(-v)] for k,v in sorted(row.items()) if k>=count]}
                            for mon,row in sorted(original.items())]}
(HERE/(prefix+'_product_map.json')).write_text(json.dumps(product,indent=2)+'\n')
receipt={'exact_all_coefficient_equations_pass':True,'matching_rank':len(pivots),'coefficient_equations':len(original),
         'raw_gram_variables':count,'auxiliaries':len(free),'blocks':len(sizes),'order':offset,'max_block_order':max(sizes),
         'rational_upper_records':len(entries),'full_set_adequacy_asserted':False}
(HERE/(prefix+'_pencil_verification.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
