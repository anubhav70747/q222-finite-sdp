"""Independent bounded exact full-Z6 reflection support check; no solver."""
import hashlib,json,pathlib,sys
import sympy as sp
R,S=sp.symbols('R S',nonzero=True)
W,X,Y,Z=sp.symbols('W X Y Z')
psi=(W,X,Y,Z)
# Independent original full state congruence, with u=R^6 and v=S^6.
Gpsi=(S**6*W-X,-W-S**6*X,-S**6*Y+Z,Y+S**6*Z)
base=[(i,j,l) for i in range(7) for j in range(7) for l in range(4)]
ext=[(i,j,l) for i in range(7) for j in range(13) for l in range(4)]
lookup={triple:index for index,triple in enumerate(ext)}
# Entry tuple = (Bob exponent shift, output state coordinate, sign).
Gterms={0:((6,0,1),(0,1,-1)),1:((0,0,-1),(6,1,-1)),2:((6,2,-1),(0,3,1)),3:((0,2,1),(6,3,1))}
columns=[]
for i,j,l in base:
    block=[[lookup[(6-i,j+shift,ell)],sign] for shift,ell,sign in Gterms[l]]
    actual=sp.expand(R**6 * (R**-1)**i * S**j * Gpsi[l])
    rebuilt=sp.expand(sum(sign*R**ext[index][0]*S**ext[index][1]*psi[ext[index][2]] for index,sign in block))
    assert sp.expand(actual-rebuilt)==0,(i,j,l)
    # Independent monomial-coefficient extraction, avoiding negative-power artifacts.
    p=sp.Poly(actual,R,S,W,X,Y,Z)
    expected={}
    for index,sign in block:
        ii,jj,ll=ext[index]
        key=(ii,jj)+tuple(int(h==ll) for h in range(4))
        expected[key]=sign
    assert dict(p.terms())==expected,(i,j,l,p.terms(),expected)
    assert len(block)==2 and len({index for index,_ in block})==2
    columns.append(block)
assert len(base)==196 and len(ext)==364 and sum(map(len,columns))==392
assert all(sign in (-1,1) and 0<=index<364 for block in columns for index,sign in block)
# In column convention h6=L6.T*zext, Q>=0 transports as L6 Q L6.T>=0.
mapping={
    'description':'Full-Z6 fixed reflection transport, h6(u,v,psi)=u*z6(1/u,v,G(v)*psi)=L6.T*zext',
    'formal_variables':'R=u**(1/6), S=v**(1/6)',
    'state_order':['W','X','Y','Z'],
    'base_features':[list(x) for x in base],
    'extended_features':[list(x) for x in ext],
    'L6_shape':[364,196],
    'L6_sparse_columns':columns,
    'nonzero_entries':392,
    'coefficient_ring':'integers; every nonzero entry is +1 or -1',
    'identity_checks':196,
    'independent_polynomial_coefficient_checks':196,
    'scope':'Support/transport identity only. No PSD admission, native domain, closure or global theorem check.'
}
p=pathlib.Path(__file__).with_name('full_Z6_reflection_sparse_map.json')
p.write_text(json.dumps(mapping,indent=2,sort_keys=True)+'\n',encoding='utf8')
print('PASS all 196 formal full-Z6 reflection identities in R,S,W,X,Y,Z')
print('PASS independent exact polynomial monomial-coefficient checks for all 196 rows')
print('PASS fixed L6 matrix: 364 x 196 with exactly 392 signed integer entries')
print('PASS sparse columns map (i,j,l) to (6-i,j+shift,ell), shift in {0,6}')
print('CONDITIONAL CONSEQUENCE: any oriented source Gram in any fixed library inside span(Z6), including Z18 when so embedded, transports into K6,12 via this same fixed reflection construction.')
print('No universal C1/C2 PSD assertion checked; no solver or sweep run.')
print('JSON SHA256 '+hashlib.sha256(p.read_bytes()).hexdigest())
print('Runtime Python '+sys.version.split()[0]+'; SymPy '+sp.__version__)
