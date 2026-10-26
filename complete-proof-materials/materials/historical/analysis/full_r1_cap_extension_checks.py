from fractions import Fraction as F
import json
from pathlib import Path
d=F(1,1024); S=F(81,65536); V=F(1,256);e=2*d;D=e*e+S*S
f0=1-d-d**3;R0=1+2*d+2*d*d+4*d**3+d**4*(1+2*d)
rPsi=32*R0/f0*S+F(512,3)*(d*d+d*S+S*S)
rR=F(64,3)*R0*(1+d*d)/f0*S
rows={
'Ru':[17*rPsi,17*(64*d*S+F(256,3)*S*S),17*F(1280,3)*(e+S)**2,17*rR],
'Jalpha':[16*(8*d*e*e+F(32,3)*e*e*S),0,0,16*F(64,3)*e*e*S],
'O':[16*F(128,3)*d*V*(1+2*S)*D,16*32*V*V*D,16*F(512,3)*d*S*V*D,16*F(64,3)*d*V*D],
'kc':[12*d+48*d*S,18*V,96*d*S,6*d],
'H1_H2':[256*(d*d*V+V*V+d**3*S+d*S*V),0,256*(V*V+d*S*V),0],
'gc':[2*(d*d+V),0,0,0]}
tot=[sum(a[i] for a in rows.values()) for i in range(4)]
lim=[F(3,4),F(1,8),F(1,8),F(1,2)]
assert all(a<b for a,b in zip(tot,lim))
core=F(1,16)-d-F(9,4096)-V-S*V
assert core>F(1,32)
out={'source_delta':str(d),'cap_u':'9/256','cap_v':'1/16','rPsi_bound':str(rPsi),'rR_bound':str(rR),'loss_totals':list(map(str,tot)),'loss_total_decimals':list(map(float,tot)),'strict_upper_bounds':list(map(str,lim)),'comparisons':[a<b for a,b in zip(tot,lim)],'core_reserve':str(core),'core_reserve_gt_1_32':core>F(1,32),'rows':{k:list(map(str,v)) for k,v in rows.items()}}
p=Path('analysis/full_r1_cap_extension_checks.json');p.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))


