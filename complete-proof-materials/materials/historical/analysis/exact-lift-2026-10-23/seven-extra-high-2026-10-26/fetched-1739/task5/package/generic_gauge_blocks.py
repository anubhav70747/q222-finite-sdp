import sympy as s,sys
sys.path.insert(0,'/mnt/data/Q222_C1_fixed_admission_2026_10_26')
from eight_gauge_constructor import *
from c1_small_map import *
d,m,C,A,D,al,be=s.symbols('d m C A D alpha beta')
uc,vc,wc,xc,yc,zc=s.symbols('uc vc wc xc yc zc',nonzero=True)
raw=sym_parent((d,0,0,m,C,A,D,al,be))
R=coefficients(raw.subs({u:uc*u,v:vc*v,W:wc*W,X:xc*X,Y:yc*Y,Z:zc*Z}))
r=R[INVARIANT_ROWS,:]
base=s.zeros(44,1)
for j,x in zip(PIV_COLS,FIXED_INVERSE*r[PIV_ROWS,:]):base[j]=s.factor(x)
sv=s.symbols('s0:8');Z=base+sum((sv[j]*GAUGES[j] for j in range(8)),s.zeros(44,1))
H=s.zeros(12)
for j,col in enumerate(BLKCOL):
 if col==78: continue
 p,q=PAIR[col];H[p,q]=H[q,p]=Z[j]
with open('/mnt/data/Q222_C1_ORIENTED_STRICT_2026_10_26/generic_gauge_blocks.txt','w') as file:
 for name,inds in [('H7',range(7)),('H5',range(7,12))]:
  file.write(name+'\n')
  for i in inds:
   for j in inds:
    if i>j:continue
    file.write(f'{i},{j}: '+str(s.factor(H[i,j]))+'\n')
 file.write('odd: '+str(Z[-1])+'\n')
