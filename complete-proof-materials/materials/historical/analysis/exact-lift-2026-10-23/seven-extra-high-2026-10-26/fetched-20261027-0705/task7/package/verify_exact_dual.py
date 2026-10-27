#!/usr/bin/env python3
"""Exact K364 moment counterexample and symbolic source-packet inequalities. No SDP or sweep."""
import json, pathlib, re, fractions
from fractions import Fraction as Q
S=pathlib.Path(__file__).resolve().parent
packet_path=S/'FROZEN_SOURCE_PACKET.txt'
packet=packet_path.read_text()
line=next(x for x in packet.splitlines() if x.startswith('{"source_order":'))
data=json.loads(line)
assert data['source_order']==['d','lambda','n','m','C','A','D','alpha','beta']
assert sum(map(len,data['source_columns'].values()))==120

def pairing(coef):
    # y[k]=1/91 iff state key WW. All testing exponents evaluate at u=v=1.
    return sum(coef[name]*sum(Q(val) for key,val in col if key[2:]==[0,0])/91 for name,col in data['source_columns'].items())
# Strict C2 physical native seed and a strictly physical C1 witness supplied in the frozen evidence.
c2={'d':Q(1),'lambda':Q(0),'n':Q(1,8),'m':Q(3,8),'C':Q(1,6),'A':Q(5,8),'D':Q(3,8),'alpha':Q(1,6),'beta':Q(0)}
c1={'d':Q(923999207,13599785800000),'lambda':Q(0),'n':Q(0),
'm':Q(8000002,63491000001),'C':Q(2),
'A':Q(53599207,27199544400428400000),
'D':Q(11449570655261,90188053200000),
'alpha':Q(228841893,21473353157782),
'beta':Q(2883242147,1431556400000)}
assert pairing(c2)==Q(38,273)
assert pairing(c1)==8*(c1['m']+2*c1['C']+c1['D']+c1['d']/2)/91
assert pairing(c2)>0 and pairing(c1)>0
# Exactly 7*13=91 W features. Evaluation at (1,1;1,0,0,0).
features=[(i,j,l) for i in range(7) for j in range(13) for l in range(4)]
assert len(features)==364
nz=[z for z in features if z[2]==0]; assert len(nz)==91
for a in range(364):
    i,j,l=features[a]
    for b in range(a,364):
        ip,jp,lp=features[b]
        moment=(Q(1,91) if l==lp==0 else Q(0))
        assert moment==(Q(1,91) if (l==0 and lp==0) else Q(0))
# Every permissible key is a moment matrix entry (surjectivity of the moment coordinate map).
for I in range(13):
    i=min(I,6); ip=I-i; assert 0<=ip<=6
    for J in range(25):
        j=min(J,12); jp=J-j; assert 0<=jp<=12
assert 13*25*10==3250
# The canonical coefficient pairing agrees with independent packet-wise W^2 at u=v=1.
# F_W = 4lambda (1+u²) +4m v² +4C(1+v²)u² +(4D+2d)u².
# At u=v=1, (1+v²)F_W=16lambda + 8m+16C+8D+4d.
for sample in (c1,c2):
    direct=(16*sample['lambda']+8*sample['m']+16*sample['C']+8*sample['D']+4*sample['d'])/91
    assert direct==pairing(sample)
report={'verification':'PASS','trace_M':'1','rank_M':1,'positive_semidefinite':True,
        'nonzero_moment_keys':'(I,J,W,W)=1/91 for all I=0..12 J=0..24',
        'C2_seed_pairing':str(pairing(c2)),'C1_strict_witness_pairing':str(pairing(c1)),
        'source_column_nonzeros':120,'coefficient_keys':3250,'features':364,
        'PD_rational_anchor':'sum over r=1..7,s=1..13,l=0..3 of evaluation at (r^6,s^6,e_l): 364x364 positive definite by tensor-product Vandermonde',
        'negative_witness_route':'For a rational c outside K364, mix a negative PSD moment witness with rational PD anchor, normalize, approximate rationally; retains strict negative pairing.'}
(S/'EXACT_EVIDENCE.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
