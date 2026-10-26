"""Check the fixed order-33 counterexample against the published COO arrays.

Only standard-library exact arithmetic; no compiler import or SDP solver.
The supplied LDL factors and supplied source-pairing formula are not used.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
PACKAGE = HERE.parent/'pro-order33-obstruction'/'Q222_Halfline_Counterexample'
PUBLIC = PROJECT/'github-q222-finite-sdp'/'extracted-lifts'
cert_path = PACKAGE/'data'/'dyadic_order33_separator.json'
cert_bytes = cert_path.read_bytes()
cert = json.loads(cert_bytes)
source_order = ['d','lambda','n','m','C','A','D','alpha','beta']
assert cert['source_order'] == source_order
moments = [[[list(map(Q,row)) for row in A] for A in side] for side in cert['original_moments']]
hatted = [[[list(map(Q,row)) for row in A] for A in side] for side in cert['hatted_moments']]
printed= re.findall(r'\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}',
                   (PACKAGE/'proof'/'dual_six_matrices.tex').read_text(),re.S)
assert len(printed)==6
def read_cell(cell):
    cell=cell.strip()
    power=re.fullmatch(r'(-?\d+)/2\^\{?(\d+)\}?',cell)
    return Q(int(power[1]),2**int(power[2])) if power else Q(cell)
for block,text in zip([A for side in hatted for A in side],printed):
    rows=[[read_cell(cell)/2**30 for cell in row.split('&')]
          for row in text.split(r'\\') if row.strip()]
    assert rows==block
for side in range(2):
    delta = list(map(Q,cert['state_congruences'][side]))
    for k in range(3):
        assert moments[side][k] == [[Q(16)**k*delta[i]*hatted[side][k][i][j]*delta[j]
                                    for j in range(4)] for i in range(4)]
for k in range(3):
    assert all(moments[0][k][i][j] == moments[1][k][i][j] for i in range(2) for j in range(2,4))

def block(U,R,W):
    return [U[i]+R[i] for i in range(4)]+[R[i]+W[i] for i in range(4)]

def ldl_pivots(A):
    assert A == [list(row) for row in zip(*A)]
    B = [row[:] for row in A]
    out=[]
    for k in range(len(B)):
        pivot = B[k][k]
        assert pivot > 0, (k,pivot)
        out.append(pivot)
        for i in range(k+1,len(B)):
            for j in range(k+1,len(B)):
                B[i][j] -= B[i][k]*B[k][j]/pivot
    return out

dual_blocks = [block(*moments[0]),block(*moments[1]),moments[0][1],moments[1][1]]
pivots = [ldl_pivots(A) for A in dual_blocks]
W = [[Q(0) for j in range(33)] for i in range(33)]
offset=0
for A in dual_blocks:
    for i,row in enumerate(A):
        for j,value in enumerate(row):
            W[offset+i][offset+j]=value
    offset+=len(A)
assert offset==24

ell=[-23803618,906393494300170,2744379904,1378590360912,
     5935589516,20425,2600028,91248,278104928]
expected_coeff=[Q(a,2**42) for a in ell]+[Q(0)]*44
coo_results=[]
for N in [24,33]:
    path=PUBLIC/'canonical-completion'/'halfline'/f'HALFLINE_QUADRATIC_SELECTOR_ORDER{N}.coo.json'
    raw=path.read_bytes(); data=json.loads(raw)
    assert data['source_order']==source_order and data['auxiliary_count']==44
    assert data['blocks']==[8,8,4,4]+([1]*9 if N==33 else [])
    assert data['indexing']=='zero-based; full symmetric COO; do not mirror or double entries'
    coefficients=[Q(0)]*len(data['parameter_order'])
    for k,i,j,value in data['entries']:
        coefficients[k]+=W[j][i]*Q(value)
    assert coefficients==expected_coeff
    coo_results.append({'file':str(path.relative_to(PROJECT)),
                        'sha256':hashlib.sha256(raw).hexdigest(),
                        'order':N,'coefficient_count_checked':len(coefficients),
                        'all_44_auxiliary_pairings_zero':True,
                        'source_pairings':list(map(str,coefficients[:9]))})

cstar=list(map(Q,cert['source']))
cbar=[Q(1),Q(1,2**34),Q(1,2**14),Q(1,2**24),Q(1,576),Q(161),Q(1,2),Q(76),Q(1,250)]
public_source_path=PUBLIC/'outer-candidates'/'strict_outer_physical_source.json'
public_source_raw=public_source_path.read_bytes()
public_source=json.loads(public_source_raw)
assert [cstar[i] for i in [0,3,4,5,6,7,8]]==list(map(Q,public_source['source']))
assert cstar[1:3]==[0,0]
assert cbar[0]==cstar[0] and all(a>0 and a>=b for a,b in zip(cbar,cstar))
negative=sum(a*b for a,b in zip(ell,cbar))
assert negative==Q(-5417798011632274183,9663676416000)<-560635
old_negative=sum(a*b for a,b in zip(ell,cstar))
assert old_negative==Q(-303722700569565123211387668040,267603838100324717114473)<0

# Match the printed four elementary AM-GM brackets to the independently
# expanded, six-term Laurent determinant retained by the source verifier.
validity=json.loads((HERE/'pro_completion_selfcontained_data.json').read_text())
rho1=1267917833244593541025180323231
rho2=762274793296209851400745596627083
rho3=2737791187307806620394799996019
rho4=79165812050999420348852524500
expected_weights={(-1,0):rho1+rho4,(-1,1):rho3,
                  (0,-1):rho2+rho3+rho4,(0,1):rho2,
                  (1,0):rho1+rho3,(1,1):rho4}
assert {tuple(row['exponent']):int(row['weight']) for row in validity['integer_weights']}==expected_weights
assert all(x>0 for x in [rho1,rho2,rho3,rho4])

result={'certificate_sha256':hashlib.sha256(cert_bytes).hexdigest(),
        'published_COO_checks':coo_results,
        'four_dual_blocks_strictly_positive_definite':True,
        'direct_rational_LDL_pivots':[[str(x) for x in p] for p in pivots],
        'congruence_from_all_six_hatted_matrices_exact':True,
        'all_six_printed_TeX_matrices_match_data':True,
        'all_twelve_selector_overlap_equations_exact':True,
        'four_AMGM_brackets_match_direct_determinant_weights':True,
        'source_order':source_order,'cbar':list(map(str,cbar)),
        'cbar_minus_cstar':list(map(str,[a-b for a,b in zip(cbar,cstar)])),
        'cstar_matches_published_source':True,
        'published_valid_source_sha256':hashlib.sha256(public_source_raw).hexdigest(),
        'integer_functional':ell,'integer_functional_value_at_cbar':str(negative),
        'actual_dual_pairing_at_cbar':str(negative/Q(2**42)),
        'integer_functional_value_at_cstar':str(old_negative),
        'scope':'Exact infeasibility of the canonical order24/33 half-line quadratic selector, extending to every polynomial selector by its forced degree<=2. All-state source validity is established by the separately reviewed determinant/AM-GM proof. No conclusion about a relabelled global assembly or the enlarged ramified cone is made.'}
(HERE/'order33_counterexample_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='direct_rational_LDL_pivots'},indent=2))
