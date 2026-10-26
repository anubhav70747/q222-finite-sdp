#!/usr/bin/env python3
"""Deterministic *exact* C2 fixed-Z6 support-face coefficient and dual-moment map.

Requirements: Python 3, SymPy. No solver, numeric samples, network, or external repository.
Run: python c2_fixed_map.py --check
Optional: python c2_fixed_map.py --write-json c2_map_data.json

The full native source has six strictly positive events (n,m,C,A,D,alpha),
d=1, lambda=beta=0. The fixed support face retains 154 sixth-root features.
The contact is scaled to U=V=phi_W=phi_X=phi_Y=phi_Z=1 before coefficient
comparison. Four sign-parity Gram blocks are compressed by their contact
vectors. This program does NOT assert a PSD Gram for arbitrary source inputs.
"""
import argparse, json
from collections import defaultdict
import sympy as sp

a,b,c,t,r=sp.symbols('a b c t r',positive=True)
n,m,C,A,D,al=sp.symbols('n m C A D alpha')
u,v,W,X,Y,Z=sp.symbols('u v W X Y Z')
STATES=(W,X,Y,Z)
PARITIES=((0,0),(0,1),(1,0),(1,1))

def native_weights():
    z=c**2; Ac=1+z*r
    h=b-a*a*t; g=(1+z)*b-(1-r+a*a*z*t); j=1+z*r-t
    delta=(1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r+a*a*z*((1+z)*b+1+z*r*(2-r)))
    N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
    wn=((b+r)*Ac*(Ac+a*a*z*t)-r*delta)/(2*r*delta)
    wm=(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*delta)
    wC=Ac*N/(2*a*a*t*delta)
    wA=a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*delta)
    wD=(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*delta)
    wal=z*z*(b+r)*N/(2*(1-z*h)*delta)
    return {n:wn,m:wm,C:wC,A:wA,D:wD,al:wal}, {'h':h,'g':g,'j':j,'zh':1-z*h,'bz':b-a*a*z*r,'Delta':delta,'N':N,'n_numer':sp.factor((b+r)*Ac*(Ac+a*a*z*t)-r*delta)}

def is_strict_native(parameters):
    vals={a:sp.Rational(parameters[0]), b:sp.Rational(parameters[1]),
          c:sp.Rational(parameters[2]), t:sp.Rational(parameters[3]),r:sp.Rational(parameters[4])}
    if not all(vals[s] > 0 for s in (a,b,c,t,r)) or not vals[r]<1: return False
    ww,gaps=native_weights()
    return (all(sp.cancel(q.subs(vals))>0 for q in gaps.values()) and
            all(sp.cancel(q.subs(vals))>0 for q in ww.values()))

def full_parent():
    du=1+u*u; dv=1+v*v
    return sp.expand(4*n*du*dv*X*X+4*m*(v*W-X+u*v*Y-u*Z)**2+
        4*C*dv*(u*W-Y)**2+4*A*du*dv*Y*Y+
        (4*D+2)*(u*W+u*v*X-Y-v*Z)**2+
        (4*al+2)*du*(v*Y-Z)**2+2*dv*(X+u*Z)**2-2*du*dv*Z*Z)

def target_template():
    """Full 4-state SG2 coefficient tensor, with all mixed normal terms."""
    p=sp.Poly(full_parent(),u,v,W,X,Y,Z)
    out={}
    for powers,co in p.terms():
        alpha,beta,*state=powers
        assert sum(state)==2
        pair=tuple(i for i,z in enumerate(state) for _ in range(z))
        assert len(pair)==2 and 0<=alpha<=2 and 0<=beta<=2
        key=(6*alpha,6*beta,min(pair),max(pair))
        out[key]=sp.expand(co)
    return out

def target_native():
    """Exact rational target at standard contact, no testing inverses."""
    wt,_=native_weights()
    uc=c; vc=a*c; psi=(b/a,-c*r,c*(b/a-a*t),sp.Integer(1))
    out={}
    for key,co in target_template().items():
        p,q,l,m_=key
        res=co.subs(wt,simultaneous=True)*uc**(p//6)*vc**(q//6)*psi[l]*psi[m_]
        out[key]=sp.cancel(res)
    return out

def keep(i,j,l):
    return not ((l==0 and i+j<6) or (l==3 and i+j>6))

def blocks():
    return [[(i,j,l) for l in range(4) for i in range(7) for j in range(7)
             if (i%2,j%2)==par and keep(i,j,l)] for par in PARITIES]

def product_key(v1,v2):
    i,j,l=v1; ii,jj,ll=v2
    return (i+ii,j+jj,min(l,ll),max(l,ll))

def diff_product(f1,ref1,f2,ref2):
    terms=defaultdict(int)
    for z1,s1 in ((f1,1),(ref1,-1)):
        for z2,s2 in ((f2,1),(ref2,-1)):
            terms[product_key(z1,z2)]+=s1*s2
    return {key:val for key,val in terms.items() if val}

def coefficient_map():
    """Exact integer sparse map: rows 412, sym 4-block variables 3021."""
    bl=blocks()
    row_keys=sorted({product_key(fi,fj) for block in bl for fi in block for fj in block})
    row_index={key:i for i,key in enumerate(row_keys)}
    cols=[]; columns=[]; adj_blocks=[]
    for bi,block in enumerate(bl):
        pivot=block[0]; rem=block[1:]
        adj_blocks.append((pivot,rem))
        for j in range(len(rem)):
            for i in range(j+1):
                coef=diff_product(rem[i],pivot,rem[j],pivot)
                factor=1 if i==j else 2
                col={row_index[k]:factor*v for k,v in coef.items()}
                cols.append(col); columns.append((bi,i,j))
    return {'blocks':bl,'row_keys':row_keys,'columns':columns,'sparse_columns':cols,'adj_blocks':adj_blocks}

def contact_annihilators(row_keys):
    # rows are angular first derivatives (scaled by 6) and 4 state derivatives
    N=[]
    for q in range(4): N.append([int(l==q)+int(m==q) for _,_,l,m in row_keys])
    N.extend([[p for p,_,_,_ in row_keys],[q for _,q,_,_ in row_keys]])
    return N

def verify_map():
    M=coefficient_map(); rows=M['row_keys']; columns=M['sparse_columns']
    orders=list(map(len,M['blocks']))
    assert orders==[52,36,36,30],orders
    assert len(rows)==412,len(rows)
    assert len(columns)==3021,len(columns)
    N=contact_annihilators(rows)
    assert sp.Matrix(N).rank()==6
    assert all(sum(N[k][row]*val for row,val in col.items())==0
               for k in range(6) for col in columns)
    src=target_template()
    assert all(key in rows for key in src)
    print('EXACT: parity blocks',orders,'rows',len(rows),'columns',len(columns))
    print('EXACT: all columns annihilated by six independent contact derivatives')
    print('EXACT: all full-parent coefficient keys supported by fixed Z6 face')
    print('RANK INPUT: full 412x3021 coefficient map has exact rank 406')
    print('  proof: native contact_rank.py: independent rank-406 modular minor, upper bound via six derivatives')
    return M

def verify_native_contact():
    M=coefficient_map()
    nt=contact_annihilators(M['row_keys'])
    templ=target_template(); wt,_=native_weights()
    uc=c; vc=a*c; psi=(b/a,-c*r,c*(b/a-a*t),sp.Integer(1))
    # Build six source-affine contact polynomials with compact symbolic expressions
    for k in range(6):
        value=0
        for key,pco in templ.items():
            p,q,l,m_=key
            ann=nt[k][M['row_keys'].index(key)]
            if ann:
                value += ann*pco*uc**(p//6)*vc**(q//6)*psi[l]*psi[m_]
        val=sp.cancel(value.subs(wt,simultaneous=True))
        assert val==0, (k,sp.factor(val))
    print('EXACT: native source substitution annihilated by all six full-parent contact conditions')

def write_data(path):
    M=coefficient_map(); N=contact_annihilators(M['row_keys'])
    trace=[0]*len(M['row_keys'])
    for bi,(pivot,rem) in enumerate(M['adj_blocks']):
        for feat in rem:
            for key,v in diff_product(feat,pivot,feat,pivot).items():
                trace[M['row_keys'].index(key)]+=v
    o={'schema':'C2-face-Z6-contact-reduced-v1',
       'parameter_order':['a','b','c','t','r'],
       'block_orders':[len(x) for x in M['blocks']],
       'support_blocks':M['blocks'],
       'coefficient_rows':M['row_keys'],
       'reduced_symmetric_columns':M['columns'],
       'sparse_coefficient_columns':[[[i,v] for i,v in sorted(c.items())] for c in M['sparse_columns']],
       'contact_annihilators':N,'dual_trace_vector':trace,
       'known_rank_exact':406, 'rank_evidence':'independent-review/checks/gram-lift/contact_rank.py, 65521 modular pivot; not independently replayed by this generator'}
    with open(path,'w') as f: json.dump(o,f,separators=(',',':'))
    print('EXACT: wrote sparse self-contained coefficient/contact/dual data to',path)

if __name__=='__main__':
    pa=argparse.ArgumentParser(); pa.add_argument('--check',action='store_true');pa.add_argument('--symbolic-contact',action='store_true');pa.add_argument('--write-json'); args=pa.parse_args()
    if args.check: verify_map()
    if args.symbolic_contact: verify_native_contact()
    if args.write_json: write_data(args.write_json)
    if not any((args.check,args.symbolic_contact,args.write_json)): pa.print_help()
