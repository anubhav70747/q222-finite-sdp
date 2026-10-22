#!/usr/bin/env python3
"""Independent exact audit of the newly evaluated weighted family and pencil.

This checker imports neither compiler and uses SymPy's Groebner reduction
from the original six cover equations.  It does not read the old order-669
pencil.  Python 3.10+ and SymPy suffice.  Optional --full-source-identity
checks all ten state coefficients symbolically in all nine source variables.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix


HERE = Path(__file__).resolve().parent
NAMES = (
    "WEIGHTED_FIXED_FUNCTIONS.json",
    "WEIGHTED_GRAM_PRODUCT_MAP.coo.json",
    "WEIGHTED_CONE_ORDER133.coo.json",
    "CANONICAL_SOURCE_QUARTIC_MAP.json",
)


def load(name):
    raw = (HERE/name).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def sparse(data):
    rows, cols = data["shape"]
    out = {}
    for i,j,p,q in data["entries"]:
        assert 0 <= i < rows and 0 <= j < cols and q > 0
        assert (i,j) not in out
        out[i,j] = s.Rational(p,q)
    return s.SparseMatrix(rows,cols,out)


def monomial(exponents,variables):
    return s.prod(v**e for v,e in zip(variables,exponents) if e)


def assert_zero(expr,description):
    assert expr == 0, (description,str(expr)[:1200])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--full-source-identity",action="store_true")
    args=parser.parse_args()
    start=time.time()
    values=[load(name) for name in NAMES]
    feature,product,pencil,quartic=[v[0] for v in values]
    hashes={name:value[1] for name,value in zip(NAMES,values)}
    variables=s.symbols("a0 a1 b0 b1 x0 x1 y0 y1 W X Y Z")
    a0,a1,b0,b1,x0,x1,y0,y1,W,X,Y,Z=variables
    anglevars=variables[:8]
    statevars=variables[8:]
    assert list(map(str,variables))==feature["variables"]==product["variables"]
    relations=[a0*a0+a1*a1-1,b0*b0+b1*b1-1,
               x0**6-a0*a0,x1**6-a1*a1,y0**6-b0*b0,y1**6-b1*b1]
    gb=s.groebner(relations,*variables,order="lex",domain=s.QQ)
    assert len(gb.polys)==6

    @lru_cache(maxsize=None)
    def nf(expr):
        return s.expand(gb.reduce(s.expand(expr))[1])

    # The mathematical row table is entered independently of the JSON.
    R=s.Matrix([[b0,b1,0,0],[0,1,0,0],
                [a0*b1,-a0*b0,a1*b1,-a1*b0],
                [a1,0,-a0,0],[0,0,1,0],
                [a1*b0,a1*b1,-a0*b0,-a0*b1],
                [0,0,b1,-b0],[0,a0,0,a1]])
    rowlocals=dict(zip(map(str,variables),variables))
    assert all(s.expand(s.sympify(feature["rows"][i][j],locals=rowlocals)-R[i,j])==0
               for i in range(8) for j in range(4))
    psi=s.Matrix(statevars)
    original=[]
    original_names=[]
    pairs={}
    for i,j in itertools.combinations(range(8),2):
        pairs[i,j]=[]
        for k in range(8):
            if k in (i,j):
                continue
            det=s.Matrix.vstack(R[i,:3],R[j,:3],R[k,:3]).det()
            pairs[i,j].append(len(original))
            original.append(nf(det*(R[k,:]*psi)[0]))
            original_names.append(f"f_{i+1}_{j+1}_{k+1}")
    E=[e for e in itertools.product((-1,0,1),repeat=2) if e!=(0,0)]
    N=sorted({(3*e[0],3*e[1]) for e in E}
             |set(itertools.product((-1,1),repeat=2))
             |set(itertools.product((-1,1),(-2,2)))
             |set(itertools.product((-2,2),(-1,1)))|{(0,0)})
    assert len(N)==21
    for nu in N:
        original.append(nf(4*x0**(3-nu[0])*x1**(3+nu[0])
                           *y0**(3-nu[1])*y1**(3+nu[1])*Z))
        original_names.append(f"g_{nu[0]}_{nu[1]}")
    assert len(original)==feature["original_indexed_count"]==189
    assert original_names==feature["original_feature_names"]
    basis_indices=feature["global_basis_original_indices"]
    assert len(basis_indices)==len(feature["features"])==104
    basis=[]
    for index,data in zip(basis_indices,feature["features"]):
        polynomial=s.expand(sum(s.Rational(c)*monomial(exp,variables)
                                for exp,c in data["polynomial"]))
        assert nf(polynomial)==polynomial
        assert_zero(polynomial-original[index],"stored global feature")
        assert data["name"]==original_names[index]
        assert all(sum(exp[8:])==1 for exp in s.Poly(polynomial,*variables).monoms())
        basis.append(polynomial)
    globalmap=sparse(feature["all_original_features_in_global_basis"])
    assert globalmap.shape==(104,189)
    basis_vector=s.Matrix(basis)
    represented=globalmap.T*basis_vector
    for index in range(189):
        assert_zero(s.expand(represented[index]-original[index]),"global feature compression")

    coefficient_monomials=sorted({m for f in basis for m in s.Poly(f,*variables).monoms()})
    mi={m:i for i,m in enumerate(coefficient_monomials)}
    feature_coeffs={}
    for j,f in enumerate(basis):
        for m,c in s.Poly(f,*variables).terms():
            feature_coeffs[mi[m],j]=c
    feature_rank=DomainMatrix.from_Matrix(s.SparseMatrix(
        len(mi),104,feature_coeffs)).convert_to(s.QQ).rank()
    assert feature_rank==104
    print("PASS: all 189 original functions, 104 independent global functions",flush=True)

    blockvectors=[]
    accounted_pairs=set()
    for block in feature["blocks"]:
        B=sparse(block["basis_in_global_family"])
        T=sparse(block["original_to_block_basis"])
        order=block["order"]
        assert B.shape==(order,104) and T.rows==order
        vector=[s.expand(v) for v in B*basis_vector]
        assert len(block["basis_original_indices"])==order
        for j,index in enumerate(block["basis_original_indices"]):
            assert_zero(vector[j]-original[index],"block basis feature")
        if block["name"]=="determinant":
            assert block["n_indices"]==[list(n) for n in N]
            indices=list(range(168,189))
        else:
            _,istr,jstr=block["name"].split("_")
            i,j=int(istr)-1,int(jstr)-1
            assert (i,j) in pairs and (i,j) not in accounted_pairs
            accounted_pairs.add((i,j))
            assert block["k_indices"]==[k+1 for k in range(8) if k not in (i,j)]
            indices=pairs[i,j]
        assert T.cols==len(indices)
        represented_block=T.T*s.Matrix(vector)
        for j,index in enumerate(indices):
            assert_zero(s.expand(represented_block[j]-original[index]),"block compression")
        blockvectors.append(vector)
    missing=set(pairs)-accounted_pairs
    assert all(original[index]==0 for pair in missing for index in pairs[pair])
    sizes=[len(v) for v in blockvectors]
    assert sizes==product["psd_block_sizes"]==pencil["block_sizes"]
    assert len(sizes)==27 and sum(sizes)==pencil["matrix_order"]==133
    assert sum(n*(n+1)//2 for n in sizes)==535
    print("PASS: every block map; deleted pairs are identically zero",flush=True)

    rows=product["monomial_rows"]
    assert len(rows)==261 and len({tuple(e) for e in rows})==261
    for exp in rows:
        assert len(exp)==12 and sum(exp[8:])==2
    rowmon=[monomial(e,variables) for e in rows]
    adata={}
    columns=defaultdict(list)
    for i,j,c in product["entries"]:
        assert (i,j) not in adata and 0<=i<261 and 0<=j<535
        adata[i,j]=s.Rational(c)
        columns[j].append((i,s.Rational(c)))
    A=s.SparseMatrix(261,535,adata)
    expected_coordinates=[{"block":k,"row":i,"column":j}
                          for k,n in enumerate(sizes)
                          for i in range(n) for j in range(i,n)]
    assert product["gram_coordinates"]==expected_coordinates
    for j,coordinate in enumerate(expected_coordinates):
        k,i,l=(coordinate[t] for t in ("block","row","column"))
        expected=nf((1 if i==l else 2)*blockvectors[k][i]*blockvectors[k][l])
        actual=s.expand(sum(c*rowmon[r] for r,c in columns[j]))
        assert_zero(expected-actual,f"Gram product {j}")
        if j and j%100==0:
            print(f"Checked {j}/535 Groebner product identities",flush=True)
    print("PASS: all 535 Gram products modulo the six original relations",flush=True)

    offsets=[]
    offset=0
    for size in sizes:
        offsets.append(offset)
        offset+=size
    assert offsets==pencil["block_offsets"]
    coordinate_of_position={(offsets[c["block"]]+c["row"],
                             offsets[c["block"]]+c["column"]):j
                            for j,c in enumerate(expected_coordinates)}
    mdata={}
    for v,i,j,p,q in pencil["entries"]:
        assert 0<=v<535 and q>0 and i<=j
        assert (i,j) in coordinate_of_position
        key=(coordinate_of_position[i,j],v)
        assert key not in mdata
        mdata[key]=s.Rational(p,q)
    M=s.SparseMatrix(535,535,mdata)
    reconstruction=sparse(pencil["coefficient_reconstruction"])
    assert reconstruction.shape==(261,189)
    assert pencil["visible_coefficient_variables"]==189
    assert pencil["free_gram_variables"]==346
    assert pencil["constant_matrix"]=="zero"
    selected=pencil["selected_coefficient_rows"]
    assert len(selected)==189 and len(set(selected))==189
    assert [rows[i] for i in selected]==pencil["selected_monomials"]
    assert reconstruction.extract(selected,range(189))==s.eye(189)
    expected=reconstruction.row_join(s.zeros(261,346))
    assert A*M==expected
    free=pencil["free_gram_coordinates"]
    pivot=pencil["pivot_gram_coordinates"]
    assert len(free)==346 and len(pivot)==189
    assert sorted(free+pivot)==list(range(535))
    assert M.extract(free,range(535))==s.zeros(346,189).row_join(s.eye(346))
    # Together these two identities are a square left inverse for M.
    # Thus M is invertible, rank(A)=189, and every kernel direction is present.
    print("PASS: A*L=[coefficient reconstruction,0] and all 346 free directions",flush=True)

    src=s.symbols("d lambda n m C A D alpha beta")
    d,lam,n,m,C,AA,D,alpha,beta=src
    assert list(map(str,src))==quartic["variables"]
    assert quartic["homogeneous_degree"]==4 and quartic["output_dimension"]==189
    yexpr=[s.Integer(0) for _ in range(189)]
    qkeys=set()
    for row,exp,p,q in quartic["entries"]:
        assert 0<=row<189 and len(exp)==9 and sum(exp)==4 and min(exp)>=0 and q>0
        assert (row,tuple(exp)) not in qkeys
        qkeys.add((row,tuple(exp)))
        yexpr[row]+=s.Rational(p,q)*monomial(exp,src)
    bexpr=[s.expand(v) for v in reconstruction*s.Matrix(yexpr)]
    # A polynomial in the original half angles cannot retain odd sixth-root powers.
    for exp,coef in zip(rows,bexpr):
        if exp[4] or exp[6] or exp[5]%6 or exp[7]%6:
            assert_zero(coef,"spurious algebraic-root coefficient in target")
    kappa=[4*lam,4*n,4*m,4*C,4*AA,4*D+2*d,4*alpha+2*d,4*beta+2*d]
    P=sum((kappa[i]*R[i,:].T*R[i,:] for i in range(8)),s.zeros(4))
    P[3,3]-=2*d
    source_fixtures=[tuple(int(i==j) for i in range(9)) for j in range(9)]
    source_fixtures += [(1,1,1,1,1,1,1,1,1),
                        (2,3,5,7,11,13,17,19,23),
                        (s.Rational(2,3),0,s.Rational(3,5),2,0,1,s.Rational(5,7),0,4)]
    test_fixtures=[(s.Rational(3,5),s.Rational(4,5),s.Rational(5,13),s.Rational(12,13)),
                   (-s.Rational(3,5),s.Rational(4,5),s.Rational(5,13),-s.Rational(12,13)),
                   (0,1,1,0),(1,0,0,-1),
                   (s.Rational(5,13),-s.Rational(12,13),-s.Rational(8,17),s.Rational(15,17))]
    fixtures_checked=0
    for sv in source_fixtures:
        source_sub=dict(zip(src,sv))
        bvalues=[v.subs(source_sub) for v in bexpr]
        Psource=P.subs(source_sub)
        for tv in test_fixtures:
            angle_sub=dict(zip((a0,a1,b0,b1),tv))
            PP=Psource.subs(angle_sub)
            target=PP[:3,:3].det()*PP
            actual=s.zeros(4)
            for exp,coef in zip(rows,bvalues):
                if not coef:
                    continue
                assert exp[4]==exp[6]==0 and exp[5]%6==exp[7]%6==0
                angle_value=s.prod(tv[i]**exp[i] for i in range(4))
                angle_value*=tv[1]**(exp[5]//3)*tv[3]**(exp[7]//3)
                state=[i for i in range(4) for _ in range(exp[8+i])]
                i,j=state
                value=coef*angle_value
                if i==j:
                    actual[i,i]+=value
                else:
                    actual[i,j]+=value/2
                    actual[j,i]+=value/2
            assert actual==target
            fixtures_checked+=1
    print(f"PASS: quartic source map at {fixtures_checked} exact full-parent matrix evaluations",flush=True)

    full_source_checked=False
    if args.full_source_identity:
        gbsrc=s.groebner(relations,*anglevars,*src,order="lex",domain=s.QQ)
        H=P[:3,:3]
        detH=(H[0,0]*(H[1,1]*H[2,2]-H[1,2]*H[2,1])
              -H[0,1]*(H[1,0]*H[2,2]-H[1,2]*H[2,0])
              +H[0,2]*(H[1,0]*H[2,1]-H[1,1]*H[2,0]))
        by_state=defaultdict(lambda:s.Integer(0))
        for exp,coef in zip(rows,bexpr):
            if coef:
                by_state[tuple(exp[8:])]+=coef*monomial(exp[:8],anglevars)
        for i in range(4):
            for j in range(i,4):
                state_exp=tuple(int(k==i)+int(k==j) for k in range(4))
                difference=s.expand((1 if i==j else 2)*detH*P[i,j]-by_state[state_exp])
                remainder=gbsrc.reduce(difference)[1]
                assert_zero(remainder,f"full symbolic quartic source coefficient {i},{j}")
                print(f"PASS: symbolic full-source state coefficient {i},{j}",flush=True)
        full_source_checked=True

    for name,sha in hashes.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,"input changed during check"
    result={
        "status":"all independent exact checks passed",
        "scope":"New multiplied-slack Gram cone and quartic source map; no old669 recomputation",
        "implementation":"Independent row table; SymPy Groebner from six original relations; neither compiler imported",
        "input_sha256":hashes,
        "groebner_basis":[str(p.as_expr()) for p in gb.polys],
        "original_functions_checked":189,
        "global_feature_dimension":int(feature_rank),
        "psd_blocks":27,
        "total_psd_order":133,
        "zero_pair_blocks_removed":len(missing),
        "gram_product_columns_checked":535,
        "coefficient_monomial_rows":261,
        "product_map_nonzeros":len(adata),
        "pencil_nonzero_affine_entries":len(mdata),
        "coefficient_map_rank":189,
        "free_kernel_directions":346,
        "pencil_coefficient_identity":"A*M=[reconstruction,0]; free coordinates of M equal [0,I346]",
        "quartic_source_terms":len(qkeys),
        "rational_full_parent_matrix_evaluations":fixtures_checked,
        "source_fixtures":len(source_fixtures),
        "test_fixtures":len(test_fixtures),
        "symbolic_full_source_identity_all_ten_state_coefficients":full_source_checked,
        "floating_point_operations":0,
        "elapsed_seconds":round(time.time()-start,3),
    }
    output=HERE/"evaluated_matrix_check_results.json"
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
