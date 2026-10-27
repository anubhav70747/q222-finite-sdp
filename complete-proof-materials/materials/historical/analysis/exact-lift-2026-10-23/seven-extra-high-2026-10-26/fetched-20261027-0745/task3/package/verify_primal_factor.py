#!/usr/bin/env python3
"""Exact *witness verifier*, not a witness generator and not an H2 proof.

An SOS witness is a JSON document with either
  {"theta": ["a","b","c","t","r"], "factors": [...]}
for a genuine strict native C2 contact, or
  {"weights": {"d":"...",...,"beta":"..."}, "factors": [...]}
for a direct source-vector smoke test.
Each factor is a list of [i,j,l,coeff], meaning
   coeff * u**(i/6) * v**(j/6) * psi[l],
0 <= i <= 6, 0 <= j <= 12, 0 <= l < 4.
Coefficients are exact real sympy expressions, e.g. 'sqrt(2)/3'.
An exact factor list yields P = sum_k q_k q_k^T PSD automatically.

This checks all coefficient equations, including all source-zero rows, by
sparse accumulation. It does NOT compile the 66430-column Gram matrix.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
MAP = json.loads((HERE / 'COMPACT_K364_SOURCE_MAP.json').read_text())
WEIGHT_NAMES = MAP['source_order']

def exact(x):
    if isinstance(x,(int,s.Integer)): return s.Integer(x)
    return s.sympify(str(x), locals={'sqrt':s.sqrt}, rational=True)

def source_native(theta):
    a,b,c,t,r = [exact(x) for x in theta]
    z=c*c; Ac=1+z*r; h=b-a*a*t
    g=(1+z)*b-(1-r+a*a*z*t); j=1+z*r-t
    Delta=(1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r
           +a*a*z*((1+z)*b+1+z*r*(2-r)))
    N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
    Nn=(b+r)*Ac*(Ac+a*a*z*t)-r*Delta
    gaps={'a':a,'b':b,'c':c,'t':t,'r':r,'1-r':1-r,
          'h':h,'g':g,'j':j,'1-c2*h':1-z*h,
          'b-a2*c2*r':b-a*a*z*r,'Nstar':N,'Nn':Nn,'Delta':Delta}
    for k,v in gaps.items():
        if s.ask(s.Q.positive(v)) is not True and v.is_positive is not True:
            raise ValueError('Strict sign not proved exactly for native gap '+k)
    src={'d':s.Integer(1),'lambda':s.Integer(0),'beta':s.Integer(0),
         'n':Nn/(2*r*Delta),
         'm':(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*Delta),
         'C':Ac*N/(2*a*a*t*Delta),
         'A':a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*Delta),
         'D':(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*Delta),
         'alpha':z*z*(b+r)*N/(2*(1-z*h)*Delta)}
    contact=(c,a*c,b/a,-c*r,c*h/a,s.Integer(1))
    return {k:s.factor(v) for k,v in src.items()},contact

def polynomial_coefficients(weights):
    target=defaultdict(lambda:s.Integer(0))
    for source_name, terms in MAP['source_columns'].items():
        w=weights.get(source_name,s.Integer(0))
        if w==0:continue
        for key,co in terms:
            target[tuple(key)]+=w*exact(co)
    return {k:s.factor(v) for k,v in target.items() if s.simplify(v)!=0}

def factor_coefficients(factors):
    coeff=defaultdict(lambda:s.Integer(0))
    for row in factors:
        vec={}
        for i,j,l,co in row:
            if not (isinstance(i,int) and 0<=i<=6 and isinstance(j,int) and 0<=j<=12 and isinstance(l,int) and 0<=l<4):
                raise ValueError('A factor uses a feature outside Z_(6,12): '+repr((i,j,l)))
            key=(i,j,l)
            vec[key]=s.simplify(vec.get(key,s.Integer(0))+exact(co))
        pairs=sorted((k,v) for k,v in vec.items() if v!=0)
        for ix,(z,q) in enumerate(pairs):
            for z2,q2 in pairs[ix:]:
                key=(z[0]+z2[0],z[1]+z2[1],min(z[2],z2[2]),max(z[2],z2[2]))
                coeff[key]+=q*q2*(1 if z==z2 else 2)
    return {k:s.factor(v) for k,v in coeff.items() if s.simplify(v)!=0}

def verify(data, check_contact=True):
    if ('theta' in data)==('weights' in data):
        raise ValueError('Specify precisely one of theta or weights')
    if 'theta' in data:
        weights,contact=source_native(data['theta'])
        uc,vc,W,X,Y,Z=contact
        assert weights['lambda']==weights['beta']==0
        if check_contact:
            for ix,row in enumerate(data['factors']):
                at=0
                for i,j,l,co in row:
                    at+=exact(co)*uc**s.Rational(i,6)*vc**s.Rational(j,6)*contact[2+l]
                if s.simplify(at)!=0:
                    raise AssertionError(f'factor {ix} does not annihilate the native contact: {s.simplify(at)}')
    else:
        contact=None
        weights={k:exact(data['weights'].get(k,'0')) for k in WEIGHT_NAMES}
    target=polynomial_coefficients(weights)
    produced=factor_coefficients(data['factors'])
    mismatches=[]
    for k in sorted(set(target)|set(produced)):
        residue=s.simplify(target.get(k,s.Integer(0))-produced.get(k,s.Integer(0)))
        if residue!=0:mismatches.append((k,str(residue)))
    if mismatches:
        raise AssertionError('FAIL: coefficient mismatch including zero keys: '+repr(mismatches[:15])+'; total '+str(len(mismatches)))
    return {'status':'EXACT SOS WITNESS VERIFIED (one source only)',
            'factors':len(data['factors']),'nonzero_target_keys':len(target),
            'nonzero_produced_keys':len(produced),
            'native_contact_checked':contact is not None and check_contact,
            'source_weights':{k:str(v) for k,v in weights.items()},
            'disclaimer':'A single witness, even if valid, does not prove universal H2.'}

def selftest():
    # Check the exact coefficient convention against the elementary n-event.
    # This is a harness smoke test, NOT evidence for any strict C2 source.
    factors=[]
    for i in (0,6):
        for j,co in ((0,'2'),(6,'2*sqrt(2)'),(12,'2')):
            factors.append([[i,j,1,co]])
    out=verify({'weights':{'n':'1'},'factors':factors})
    assert out['factors']==6 and out['nonzero_target_keys']==6
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('witness',nargs='?',help='JSON factor-list certificate at one exact source')
    p.add_argument('--selftest',action='store_true',help='exercise coefficient conventions on the n event, not a native C2 claim')
    args=p.parse_args()
    if args.selftest: result=selftest()
    elif args.witness: result=verify(json.loads(Path(args.witness).read_text()))
    else: p.error('Supply a witness JSON or --selftest. This code does NOT generate an H2 Gram.')
    print(json.dumps(result,indent=2))
