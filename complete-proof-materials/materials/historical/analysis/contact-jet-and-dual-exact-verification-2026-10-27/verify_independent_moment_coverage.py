"""Small independent moment coverage check; no supplied-code imports."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json,time,math
import sympy as s
start=time.perf_counter();out=Path(__file__).resolve().parent
features=[(i,j,l) for i in range(7) for j in range(13) for l in range(4)];index={x:i for i,x in enumerate(features)}
q=[int(l==0) for i,j,l in features]
assert len(features)==364 and sum(q)==91
trace=sum(Q(x*x,91) for x in q);assert trace==1
upper=0
for a,b in combinations_with_replacement(range(364),2):
    i,j,l=features[a];ip,jp,lp=features[b]
    key=(i+ip,j+jp,min(l,lp),max(l,lp))
    functional=Q(1,91) if key[2:]==(0,0) else Q(0)
    outer_entry=Q(q[a]*q[b],91)
    assert functional==outer_entry
    upper+=1
# q is nonzero. M=qq^T/91 has exact rank one, PSD, and trace one.
assert q[index[(0,0,0)]]==1
reachable=0
for I in range(13):
 for J in range(25):
  i,ip=min(I,6),I-min(I,6);j,jp=min(J,12),J-min(J,12)
  for a,b in combinations_with_replacement(range(4),2):
   f=features[index[(i,j,a)]];g=features[index[(ip,jp,b)]]
   assert (f[0]+g[0],f[1]+g[1],min(f[2],g[2]),max(f[2],g[2]))==(I,J,a,b)
   reachable+=1
assert reachable==3250
V7=s.Matrix([[r**i for i in range(7)] for r in range(1,8)])
V13=s.Matrix([[r**i for i in range(13)] for r in range(1,14)])
d7=int(V7.det());d13=int(V13.det())
prod=lambda n:math.prod(b-a for a in range(1,n+1) for b in range(a+1,n+1))
assert d7==prod(7)>0 and d13==prod(13)>0
# The 364 evaluation matrix is V7 tensor V13 tensor I4. Its determinant is
# d7^(13*4) d13^(7*4), nonzero, so its Gram is strictly PD with integer data.
report={'result':'INDEPENDENT MOMENT/ANCHOR COVERAGE PASSED','upper_entries_matched_to_outer_product':upper,'rank_one_factor_nonzero_coordinates':sum(q),'trace':'1','rank':1,'PSD_factorization':'qq^T/91 with q nonzero','actual_moment_coordinates_reached':reachable,'Vandermonde7_determinant':str(d7),'Vandermonde13_determinant':str(d13),'anchor_evaluation_matrix_order':364,'anchor_invertibility':'det(V7 tensor V13 tensor I4)=det(V7)^52 det(V13)^28 >0','elapsed_seconds':round(time.perf_counter()-start,6),'universal_source_pairing_nonnegativity_checked':False,'negative_dual_separator_constructed':False}
(out/'independent_moment_coverage_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
