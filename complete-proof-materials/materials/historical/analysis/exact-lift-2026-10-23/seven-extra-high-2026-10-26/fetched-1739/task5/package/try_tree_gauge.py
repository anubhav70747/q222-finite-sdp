import sys,sympy as S
sys.path.insert(0,'/mnt/data/Q222_C1_fixed_admission_2026_10_26')
from eight_gauge_constructor import *
from c1_small_map import *
d,m,C,A,D,al,be=S.symbols('d m C A D alpha beta')
uc,vc,wc,xc,yc,zc=S.symbols('uc vc wc xc yc zc',nonzero=True)
R=coefficients(sym_parent((d,0,0,m,C,A,D,al,be)).subs({u:uc*u,v:vc*v,W:wc*W,X:xc*X,Y:yc*Y,Z:zc*Z}))
r=R[INVARIANT_ROWS,:]
zero=S.zeros(44,1)
for j,vv in zip(PIV_COLS,FIXED_INVERSE*r[PIV_ROWS,:]):zero[j]=vv
sv=S.symbols('s0:7');delta=S.symbols('delta');gauge=zero+sum((sv[j]*GAUGES[j] for j in range(7)),S.zeros(44,1))
rt=coefficients(u*v*(X-Y)**2)[INVARIANT_ROWS,:]
for j,vv in zip(PIV_COLS,-FIXED_INVERSE*rt[PIV_ROWS,:]):gauge[j]+=delta*vv
H=S.zeros(12)
for j,col in enumerate(BLKCOL):
 if col==78: continue
 p,q=PAIR[col];H[p,q]=H[q,p]=gauge[j]
opt=[(0,5),(1,3),(1,5),(1,6),(2,3),(2,5),(3,4),(4,5),(4,6),(5,6)]
# 2,3 is new δ, so 10 entries 8 vars
expr=[H[i,j] for i,j in opt]
sol=S.linsolve(expr,(*sv,delta));print(sol)
if sol!=S.EmptySet:
 print('solution len:',len(sol.args))
M,b=S.linear_eq_to_matrix(expr,(*sv,delta));print('rank matrix',M.rank(),'aug rank?',M.row_join(b).rank())
null=M.T.nullspace();print('left rel dims',len(null));
for t0 in null:
 print('relation:',S.factor((t0.T*b)[0]))
