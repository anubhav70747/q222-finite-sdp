"""Independent exact PD/adjoint check against the published order-49 pencil.

Standard library only. Does not import Pro's verifier/compiler, use supplied
LDL factors, or run an optimizer. Integer moments are the input certificate.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import re

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
PACKAGE=HERE.parent/'pro-ramified-obstruction'/'Q222_Ramified_Quartic_Obstruction'
PUBLIC=PROJECT/'github-q222-finite-sdp'/'extracted-lifts'
certificate_path=PACKAGE/'data'/'certificate.json'
raw=certificate_path.read_bytes(); cert=json.loads(raw)
source_order=['d','lambda','n','m','C','A','D','alpha','beta']
assert cert['source_order']==source_order
den=int(cert['denominator']); assert den==15*2**60
moments=[[[[int(x) for x in row] for row in A] for A in side]
         for side in cert['moment_numerators']]
assert len(moments)==2 and all(len(side)==5 for side in moments)
for side in moments:
    for A in side:
        assert len(A)==4 and all(len(row)==4 for row in A)
        assert A==[list(row) for row in zip(*A)]
for k in range(5):
    assert all(moments[0][k][i][j]==moments[1][k][i][j]
               for i in range(2) for j in range(2,4))

# Match the ten literal printed matrices to the machine-readable moments.
for side,filename in enumerate(['moments_g.tex','moments_h.tex']):
    text=(PACKAGE/'proof'/filename).read_text()
    printed=re.findall(r'\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}',text,re.S)
    assert len(printed)==5
    for k,body in enumerate(printed):
        parsed=[[int(cell.strip()) for cell in row.split('&')]
                for row in body.split(r'\\') if row.strip()]
        assert parsed==moments[side][k]

def hankel(side,n,shift):
    return [[Q(moments[side][i//4+j//4+shift][i%4][j%4])
             for j in range(4*n)] for i in range(4*n)]

blocks=[hankel(0,3,0),hankel(0,2,1),hankel(1,3,0),hankel(1,2,1)]

def rational_pivots(A):
    """Fresh rational Schur elimination, not the supplied Bareiss routine."""
    assert A==[list(row) for row in zip(*A)]
    M=[row[:] for row in A]; pivots=[]
    for k in range(len(M)):
        pivot=M[k][k]
        assert pivot>0,(k,pivot)
        pivots.append(pivot)
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):
                M[i][j]-=M[i][k]*M[k][j]/pivot
    return pivots

pivots=[rational_pivots(A) for A in blocks]
floors=[min(p.numerator//p.denominator for p in row) for row in pivots]
assert floors==[794,12397,67408,1527738]

pencil_path=PUBLIC/'canonical-completion'/'ramified-completion'/'data'/'ramified_order49_pencil.json'
pencil_raw=pencil_path.read_bytes(); pencil=json.loads(pencil_raw)
assert pencil['source_names']==source_order
assert pencil['order']==49 and pencil['auxiliary_dimension']==148
assert pencil['matrix_convention']=='Entries list upper-triangular [row,column,rational] with zero-based indices; reflect off-diagonal entries.'
expected_layout=[(name,i,1) for i,name in enumerate(source_order)]
expected_layout += [('A_G',9,12),('B_G',21,8),('A_H',29,12),('B_H',41,8)]
assert [(b['name'],b['offset'],b['order']) for b in pencil['blocks']]==expected_layout
assert (PACKAGE/'reused'/'ramified_order49_pencil.json').read_bytes()==pencil_raw
W=[[Q(0) for j in range(49)] for i in range(49)]
for off,A in zip([9,21,29,41],blocks):
    for i,row in enumerate(A):
        for j,x in enumerate(row):W[off+i][off+j]=x
assert all(W[i][j]==0 for i in range(9) for j in range(49))
assert cert['source_sign_dual_blocks']==['0']*9

def pair(entries):
    return sum(Q(value)*W[j][i]*(1 if i==j else 2) for i,j,value in entries)

# Check source-sign rows explicitly, not just the Gram part of the pencil.
for k,A in enumerate(pencil['source_matrices']):
    assert [(i,j,Q(value)) for i,j,value in A if i<9 or j<9]==[(k,k,Q(1))]
for A in pencil['auxiliary_matrices']:
    assert all(i>=9 and j>=9 for i,j,value in A)
aux_pairings=[pair(A) for A in pencil['auxiliary_matrices']]
assert len(aux_pairings)==148 and aux_pairings==[Q(0)]*148
ell=[-97633735551856,24154193783510329,85793705264,
     1787034292680750290,3303160925893294,16956553777510394,
     13682527570574,670221724556174,366431361010]
source_pairings=[pair(A) for A in pencil['source_matrices']]
assert source_pairings==[12*Q(a) for a in ell]
assert list(map(int,cert['integer_source_functional']))==ell
assert [a/den for a in source_pairings]==list(map(Q,cert['source_functional']))

pi=[0,4,5,3,1,2,6,8,7]
cbar=[Q(1),Q(1,2**34),Q(1,2**14),Q(1,2**24),Q(1,576),Q(161),Q(1,2),Q(76),Q(1,250)]
dagger=[cbar[i] for i in pi]
assert dagger==list(map(Q,cert['source']))
assert cbar==list(map(Q,cert['established_cbar']))
assert all(x>0 for x in dagger)
negative=sum(a*b for a,b in zip(ell,dagger))
assert negative==Q(-32607070960447356138172589,9663676416000)<0
assert negative==Q(cert['integer_functional_value'])
assert Q(cert['functional_scale'])==Q(12,den)
assert 12*negative/den==Q(cert['negative_value'])
public_star_path=PUBLIC/'outer-candidates'/'strict_outer_physical_source.json'
star=json.loads(public_star_path.read_text())
s=list(map(Q,star['source'])); star9=[s[0],Q(0),Q(0)]+s[1:]
boundary=[star9[i] for i in pi]
assert all(a>=b>=0 for a,b in zip(dagger,boundary))
boundary_negative=sum(a*b for a,b in zip(ell,boundary))
assert boundary_negative==Q(-4844338143654932884745484309348314410,802811514300974151343419)<0

result={'certificate_sha256':hashlib.sha256(raw).hexdigest(),
        'published_pencil':str(pencil_path.relative_to(PROJECT)),
        'published_pencil_sha256':hashlib.sha256(pencil_raw).hexdigest(),
        'published_and_reused_pencil_bytes_identical':True,
        'all_ten_printed_integer_moments_match_certificate':True,
        'twenty_selector_overlap_equations':True,
        'four_block_orders':[12,8,12,8],
        'four_dual_blocks_positive_definite':True,
        'fresh_rational_LDL_pivots':[[str(x) for x in row] for row in pivots],
        'minimum_integer_LDL_pivot_floors':floors,
        'all_148_auxiliary_contractions_zero':True,
        'all_nine_source_contractions_equal_12ell':True,
        'source_contractions':[str(x) for x in source_pairings],
        'all_nine_source_sign_blocks_checked':True,
        'source_sign_dual_weights':['0']*9,
        'integer_functional':ell,'source_order':source_order,
        'cdagger':[str(x) for x in dagger],
        'cdagger_strictly_positive':True,
        'ell_cdagger':str(negative),
        'actual_integer_dual_pairing_cdagger':str(12*negative),
        'actual_normalized_dual_pairing_cdagger':str(12*negative/den),
        'transported_boundary_source':[str(x) for x in boundary],
        'ell_boundary':str(boundary_negative),
        'actual_integer_dual_pairing_boundary':str(12*boundary_negative),
        'scope':'Exact exclusion of cdagger and the transported boundary source from the published ramified quartic selector cone. Analytically, PSD forces every polynomial selector in s=sqrt(v)/4 to have degree<=4, so all polynomial degrees in this same ramified coordinate are excluded. Prior validity of cbar/cstar is an input, not re-audited here. Rational selectors, finer algebraic ramifications, and larger testing families are not excluded.'}
(HERE/'ramified_obstruction_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='fresh_rational_LDL_pivots'},indent=2))
