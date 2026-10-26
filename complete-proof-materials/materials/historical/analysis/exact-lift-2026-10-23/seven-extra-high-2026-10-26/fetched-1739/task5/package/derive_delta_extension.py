import sys,sympy as S
sys.path.insert(0,'/mnt/data/Q222_C1_fixed_admission_2026_10_26')
from eight_gauge_constructor import *
from c1_small_map import *
A42=M[INVARIANT_ROWS,BLKCOL]
raw=coefficients(u*v*(X-Y)**2)[INVARIANT_ROWS,:]
sol=S.zeros(44,1)
for p,z in zip(PIV_COLS,-FIXED_INVERSE*raw[PIV_ROWS,:]):sol[p]=z
assert A42*sol==-raw
print('New delta gauge, nonzero integer entries:')
for j,val in enumerate(sol):
 if val:
  print(('odd',) if BLKCOL[j]==78 else PAIR[BLKCOL[j]],val)
print('delta target', [(KEYS[INVARIANT_ROWS[i]],val) for i,val in enumerate(raw) if val])
