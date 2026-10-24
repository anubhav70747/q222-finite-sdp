import itertools, sympy as sp
E=[(i,j) for i in (-1,0,1) for j in (-1,0,1) if (i,j)!=(0,0)]
circ=[]
for k in (2,3):
    for S in itertools.combinations(E,k):
        M=sp.Matrix([[e[0] for e in S],[e[1] for e in S]])
        ns=M.nullspace()
        if len(ns)!=1: continue
        w=ns[0]; w=w/ sp.gcd(list(w)) if all(x!=0 for x in w) else w
        if all(x>0 for x in w) or all(x<0 for x in w):
            w=sp.Matrix([abs(x) for x in w]); g=sp.gcd(list(w)); w=w/g
            circ.append((S,tuple(w)))
for c in circ: print(c)
print(len(circ))
# AM-GM identities
a,b,c,dd=sp.symbols('a b c d')
print(sp.expand(a**3+b**3+c**3-3*a*b*c-sp.Rational(1,2)*(a+b+c)*((a-b)**2+(b-c)**2+(c-a)**2)))
print(sp.expand(a**4+b**4+c**4+dd**4-4*a*b*c*dd-((a**2-b**2)**2+(c**2-dd**2)**2+2*(a*b-c*dd)**2)))
