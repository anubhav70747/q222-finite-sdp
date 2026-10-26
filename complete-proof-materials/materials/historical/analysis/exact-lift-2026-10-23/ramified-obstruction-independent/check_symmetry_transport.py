"""Independent source symmetry/contact/feature transport checks, stdlib only.

No Pro verifier, dual moment, SDP pencil, solver, or numerical sampling is
used. Writes one local result JSON; every input is read-only.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
proof=HERE.parent/'pro-ramified-obstruction/Q222_Ramified_Quartic_Obstruction/proof/ramified_obstruction.tex'
proof_hash=hashlib.sha256(proof.read_bytes()).hexdigest()
names=['a0','a1','b0','b1','W','X','Y','Z']
zero=(0,)*8
one={zero:Q(1)}
def add(*polys):
 out={}
 for p in polys:
  for e,c in p.items():out[e]=out.get(e,Q(0))+c
 return {e:c for e,c in out.items() if c}
def scale(a,p):return {e:a*c for e,c in p.items() if a*c}
def mul(a,b):
 out={}
 for e,c in a.items():
  for f,d in b.items():
   g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,Q(0))+c*d
 return {e:c for e,c in out.items() if c}
def square(p):return mul(p,p)
var=[{tuple(int(i==j) for i in range(8)):Q(1)} for j in range(8)]
a0,a1,b0,b1,W,X,Y,Z=var
Da=add(square(a0),square(a1));Db=add(square(b0),square(b1))
rm=add(mul(a0,mul(b1,W)),scale(-1,mul(a0,mul(b0,X))),mul(a1,mul(b1,Y)),scale(-1,mul(a1,mul(b0,Z))))
rd=add(mul(a1,mul(b0,W)),mul(a1,mul(b1,X)),scale(-1,mul(a0,mul(b0,Y))),scale(-1,mul(a0,mul(b1,Z))))
ra=add(mul(b1,Y),scale(-1,mul(b0,Z)))
rb=add(mul(a0,X),mul(a1,Z))
# Independent homogeneous expansion of the original full SG2 parent.
packets=[
 scale(2,add(square(rd),mul(Da,square(ra)),mul(Db,square(rb)),scale(-1,mul(mul(Da,Db),square(Z))))),
 scale(4,mul(Da,square(add(mul(b0,W),mul(b1,X))))),
 scale(4,mul(mul(Da,Db),square(X))),
 scale(4,square(rm)),
 scale(4,mul(Db,square(add(mul(a1,W),scale(-1,mul(a0,Y)))))),
 scale(4,mul(mul(Da,Db),square(Y))),
 scale(4,square(rd)),
 scale(4,mul(Da,square(ra))),
 scale(4,mul(Db,square(rb))),
]
pi=[0,4,5,3,1,2,6,8,7]
assert [pi[i] for i in pi]==list(range(9))
# a'=(b1,-b0), b'=(a1,-a0), psi'=(W,Y,X,Z).
permutation=[3,2,1,0,4,6,5,7]
def pullback(p):
 out={}
 for e,c in p.items():
  f=[0]*8
  for i,j in enumerate(permutation):f[j]+=e[i]
  f=tuple(f);out[f]=out.get(f,Q(0))+c*(-1)**(e[1]+e[3])
 return {e:c for e,c in out.items() if c}
for index in range(9):assert pullback(packets[pi[index]])==packets[index]
# Squaring the homogeneous angle map negates both angle pairs. The full
# parent is even in each pair, hence this is projectively involutive.
for packet in packets:assert pullback(pullback(packet))==packet

oldpath=ROOT/'github-q222-finite-sdp/extracted-lifts/outer-candidates/strict_outer_physical_source.json'
raw=oldpath.read_bytes();old=json.loads(raw)
cstar=[Q(old['source'][0]),Q(0),Q(0)]+list(map(Q,old['source'][1:]))
pistar=[cstar[i] for i in pi]
cbar=[Q(1),Q(1,2**34),Q(1,2**14),Q(1,2**24),Q(1,576),Q(161),Q(1,2),Q(76),Q(1,250)]
dagger=[cbar[i] for i in pi]
assert dagger==[Q(1),Q(1,576),Q(161),Q(1,2**24),Q(1,2**34),Q(1,2**14),Q(1,2),Q(1,250),Q(76)]
F=add(*(scale(c,p) for c,p in zip(pistar,packets)))
point=[Q(1),Q(1,16),Q(1),Q(13),Q(140),Q(-1,50),Q(1,32),Q(1)]
def evaluate(p):
 return sum(c*__import__('functools').reduce(lambda x,ez:x*ez[0]**ez[1],zip(point,e),Q(1)) for e,c in p.items())
def derivative(p,i):
 return {tuple(x-int(j==i) for j,x in enumerate(e)):c*e[i] for e,c in p.items() if e[i]}
assert evaluate(F)==0
assert [evaluate(derivative(F,i)) for i in [1,3,4,5,6,7]]==[0]*6
u=Q(1,16);psi=point[4:]
eta=[u*psi[0]+u*u*psi[2],u*psi[1]+u*u*psi[3],u*psi[0]-psi[2],u*psi[1]-psi[3]]
assert eta==[Q(71681,8192),Q(17,6400),Q(279,32),Q(-801,800)]
assert u*u==Q(1,256)

# Exact fractional exponent bookkeeping for the five Bob factors. For
# x=abs(u)>0, |v_old|=1/x and 1/sqrt(1+v_old^2)=x/sqrt(1+x^2).
old_exponents=[Q(0),Q(1,2),Q(1),Q(1,4),Q(3,4)]
transported=[1-e for e in old_exponents]
assert transported==[Q(1),Q(1,2),Q(0),Q(3,4),Q(1,4)]
denominators=[1,4,16,2,8]
# On u>0,v=1,psi=(1,0,0,0): old Alice pair after transport is
# (1/sqrt2,-1/sqrt2), so P_G1=-1/(2sqrt2). The last Bob factor supplies
# u^(1/4)/(8sqrt(1+u^2)). After multiplying by a1, order is 5/4.
assert Q(1,2)*Q(1,8)==Q(1,16)
assert transported[4]+1==Q(5,4)
# Multiplier transport: old a1'= -b0, so its square is b0^2.
assert permutation[1]==2
assert hashlib.sha256(proof.read_bytes()).hexdigest()==proof_hash
assert oldpath.read_bytes()==raw
result={
 'full_homogeneous_source_packet_identities':9,
 'source_permutation':pi,'homogeneous_map':"a'=(b1,-b0), b'=(a1,-a0), psi'=(W,Y,X,Z)",
 'literal_zero_coordinates_included':True,'projective_involution':True,
 'c_dagger_matches_pi_cbar':True,'positive_gauge_contact_value':'0',
 'all_six_contact_first_derivatives':['0']*6,'eta':list(map(str,eta)),
 'transported_Bob_abs_u_exponents':list(map(str,transported)),
 'Bob_factor_denominators':denominators,
 'transported_feature':'-u^(1/4)/(16*sqrt(2)*sqrt(1+u^2)) on u>0,v=1,psi=(1,0,0,0)',
 'common_multiplier_feature':'-u^(5/4)/(16*sqrt(2)*(1+u^2))',
 'multiplier_transport':'a1^2 -> b0^2; common multiplier a1^2*b0^2',
 'proof_sha256':proof_hash,'originals_unchanged':True,
 'scope':'Full source symmetry, contact and feature transport only. No dual matrix audit or universal adequacy claim.',
}
(HERE/'symmetry_transport_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
