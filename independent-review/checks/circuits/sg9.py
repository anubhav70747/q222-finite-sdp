import itertools
cyc=['lam','D','C','A','al','m','be','n']
# q_e positivity: from SG6 (t_i>0 always since d>0)
Q={(-1,-1):[('lam','al')],(-1,0):[('lam','D'),('al','m')],(-1,1):[('m','D')],
   (0,-1):[('al','A'),('lam','n')],(0,1):[('be','m'),('C','D')],(1,-1):[('n','A')],
   (1,0):[('be','n'),('C','A')],(1,1):[('be','C')]}
# careful: SG6 index q_ij -> exponent (i-1,j-1)
def active(S):
    return [e for e,ps in Q.items() if any(a in S and b in S for a,b in ps)]
def coercive(act):
    # 0 in interior of conv hull: for every direction among a fine set, some active e has positive dot
    import math
    if not act: return False
    for k in range(3600):
        th=2*math.pi*k/3600; c=(math.cos(th),math.sin(th))
        if max(e[0]*c[0]+e[1]*c[1] for e in act)<=1e-12: return False
    return True
def dist(a,b):
    i,j=cyc.index(a),cyc.index(b); d=abs(i-j); return min(d,8-d)
bad=0
for r in range(9):
    for S in itertools.combinations(cyc,r):
        S=set(S); om=[x for x in cyc if x not in S]
        compat=all(dist(a,b) in (1,2) for a,b in itertools.combinations(om,2)) and len(om)<=3
        # note: SG9 says coercive <=> every two omitted compatible
        pair_compat=all(dist(a,b) in (1,2) for a,b in itertools.combinations(om,2))
        c=coercive(active(S))
        if c!=pair_compat: bad+=1; print("mismatch",sorted(S),om,c,pair_compat)
print("mismatches",bad)
