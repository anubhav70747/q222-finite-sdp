#!/usr/bin/env python3
"""Evaluate a fixed Gram product map for the adjugate-multiplied parent.

This does NOT compile an affine lift of Q222.  The represented polynomial is
det(H(c,t))*F(c,t), of source degree four.  No old order-669 calculation is run.
All arithmetic is rational/integer; the normal form uses the displayed
sixth-root cover, not numerical interpolation.
"""
from __future__ import annotations

import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
VARNAMES = "a0 a1 b0 b1 x0 x1 y0 y1 W X Y Z".split()
VARS = sp.symbols(" ".join(VARNAMES))
a0, a1, b0, b1, x0, x1, y0, y1, W, X, Y, Z = VARS
PSI = sp.Matrix([W, X, Y, Z])
ROWS = [
    [b0, b1, 0, 0],
    [0, 1, 0, 0],
    [a0*b1, -a0*b0, a1*b1, -a1*b0],
    [a1, 0, -a0, 0],
    [0, 0, 1, 0],
    [a1*b0, a1*b1, -a0*b0, -a0*b1],
    [0, 0, b1, -b0],
    [0, a0, 0, a1],
]


def polynomial(expr):
    return {tuple(p): int(c) for p, c in sp.Poly(sp.expand(expr), *VARS).terms() if c}


def add_term(out, p, c):
    out[p] += c
    if not out[p]:
        del out[p]


def normal_monomial(p):
    """Normal form for a0²=1-x1⁶,a1²=x1⁶,x0⁶=1-x1⁶ and Bob.

    Leading monomials a0²,a1²,b0²,b1²,x0⁶,y0⁶ are coprime,
    so these six monic relations are a Groebner basis.  Every product is
    homogeneous quadratic in the state; its zero on S³ is equivalent to its
    zero on all states, hence the state-sphere relation need not be reduced.
    """
    p = list(p)
    aa, ab = p[0]//2 + p[4]//6, p[2]//2 + p[6]//6
    p[5] += 6*(p[1]//2)
    p[7] += 6*(p[3]//2)
    for i in (0, 1, 2, 3):
        p[i] %= 2
    p[4] %= 6
    p[6] %= 6
    out = {}
    for ka in range(aa+1):
        for kb in range(ab+1):
            q = p.copy()
            q[5] += 6*ka
            q[7] += 6*kb
            out[tuple(q)] = (-1)**(ka+kb)*math.comb(aa, ka)*math.comb(ab, kb)
    return out


def normal(poly):
    out = defaultdict(int)
    for p, c in poly.items():
        for q, d in normal_monomial(p).items():
            add_term(out, q, c*d)
    return dict(out)


def multiply(left, right):
    raw = defaultdict(int)
    for p, c in left.items():
        for q, d in right.items():
            add_term(raw, tuple(a+b for a,b in zip(p,q)), c*d)
    return normal(raw)


def encode(poly):
    return [[list(p), c] for p,c in sorted(poly.items())]


def encode_matrix(M):
    return {"shape":list(M.shape), "entries":[[i,j,int(M[i,j].p),int(M[i,j].q)]
        for i in range(M.rows) for j in range(M.cols) if M[i,j]]}


def main():
    features = []
    blocks = []
    for i,j in itertools.combinations(range(8), 2):
        block = {"name":f"adj_{i+1}_{j+1}", "features":[], "k_indices":[]}
        for k in range(8):
            if k in (i,j):
                continue
            determinant = sp.Matrix([ROWS[i][:3], ROWS[j][:3], ROWS[k][:3]]).det()
            expr = determinant*(sp.Matrix([ROWS[k]])*PSI)[0]
            poly = normal(polynomial(expr))
            index = len(features)
            features.append({"name":f"f_{i+1}_{j+1}_{k+1}",
                             "ijk":[i+1,j+1,k+1], "polynomial":poly})
            block["features"].append(index)
            block["k_indices"].append(k+1)
        blocks.append(block)

    E = [(i,j) for i in (-1,0,1) for j in (-1,0,1) if (i,j)!=(0,0)]
    N = {(3*i,3*j) for i,j in E}
    N |= set(itertools.product((-1,1),repeat=2))
    N |= set(itertools.product((-1,1),(-2,2)))
    N |= set(itertools.product((-2,2),(-1,1)))
    N.add((0,0))
    angular = {"name":"determinant", "features":[], "n_indices":[]}
    for n,m in sorted(N):
        expr = 4*x0**(3-n)*x1**(3+n)*y0**(3-m)*y1**(3+m)*Z
        index = len(features)
        features.append({"name":f"g_{n}_{m}", "n":[n,m], "polynomial":normal(polynomial(expr))})
        angular["features"].append(index)
        angular["n_indices"].append([n,m])
    blocks.append(angular)
    assert len(features)==189 and len(N)==21 and len(blocks)==29

    # Exact linear compression, derived from the identities, not a guessed frame.
    original_features=features
    original_blocks=blocks
    feature_monomials=sorted(set().union(*(set(f["polynomial"]) for f in features)))
    F=sp.SparseMatrix(len(feature_monomials),len(features),
        {(i,j):f["polynomial"][m] for i,m in enumerate(feature_monomials)
         for j,f in enumerate(features) if m in f["polynomial"]})
    _,global_pivots=F.rref()
    G=F[:,list(global_pivots)]
    _,global_rows=G.T.rref()
    Gsquare=G[list(global_rows),:]
    global_all=Gsquare.inv()*F[list(global_rows),:]
    assert G*global_all==F
    features=[original_features[j] for j in global_pivots]
    blocks=[]
    for block in original_blocks:
        inds=block["features"]
        A=F[:,inds]
        _,pivots=A.rref()
        if not pivots:
            continue
        chosen=[inds[j] for j in pivots]
        B=F[:,chosen]
        _,rows=B.T.rref()
        T=B[list(rows),:].inv()*A[list(rows),:]
        assert B*T==A
        new={k:v for k,v in block.items() if k!="features"}
        new["basis_original_indices"]=chosen
        new["order"]=len(chosen)
        new["basis_in_global_family"]=encode_matrix(global_all[:,chosen].T)
        new["original_to_block_basis"]=encode_matrix(T)
        new["source_gram_pushforward"]="T G_original T^T; for adj_i_j this equals kappa_i*kappa_j*(T*kappa_without_ij)*(T*kappa_without_ij)^T."
        blocks.append(new)

    coordinates, columns = [], []
    for ib, block in enumerate(blocks):
        inds=block["basis_original_indices"]
        for i in range(len(inds)):
            for j in range(i,len(inds)):
                coordinates.append({"block":ib, "row":i, "column":j})
                p=multiply(original_features[inds[i]]["polynomial"], original_features[inds[j]]["polynomial"])
                columns.append({q:(1 if i==j else 2)*c for q,c in p.items()})
    monomials=sorted(set().union(*(set(p) for p in columns)))
    row_index={p:i for i,p in enumerate(monomials)}
    entries=[[row_index[p],j,c] for j,col in enumerate(columns) for p,c in sorted(col.items())]
    data={
      "schema":"Q222.adjugate-weighted-fixed-family.v2",
      "status":"Fixed Gram cone for the multiplied slack det(H)*F; not an affine Q222 lift.",
      "variables":VARNAMES,
      "domain_equalities":["a0^2+a1^2=1","b0^2+b1^2=1",
          "x0^6=a0^2","x1^6=a1^2","y0^6=b0^2","y1^6=b1^2",
          "W^2+X^2+Y^2+Z^2=1"],
      "domain_inequalities":[],
      "root_sign_convention":"All real sixth roots are allowed. Each angular circuit Gram product is even under independent sign changes of x0,x1,y0,y1.",
      "source_coefficients":["4lambda","4n","4m","4C","4A","4D+2d","4alpha+2d","4beta+2d"],
      "rows":[[str(v) for v in row] for row in ROWS],
      "original_indexed_count":len(original_features),
      "global_basis_original_indices":list(global_pivots),
      "all_original_features_in_global_basis":encode_matrix(global_all),
      "original_feature_names":[f["name"] for f in original_features],
      "features":[{**{k:v for k,v in f.items() if k!="polynomial"},
                   "polynomial":encode(f["polynomial"])} for f in features],
      "blocks":blocks,
      "source_gram_rule":"For each adj_i_j use the explicitly pushed-forward rank-one source block. The determinant block follows the twelve evaluated rational balanced-circuit matrices and source-only algebraic prescription.",
    }
    product={"schema":"Q222.exact-block-Gram-product-map.v1",
       "status":data["status"], "variables":VARNAMES,
       "normal_form_rules":["a0^2=1-x1^6","a1^2=x1^6","b0^2=1-y1^6","b1^2=y1^6","x0^6=1-x1^6","y0^6=1-y1^6"],
       "state_rule":"All terms are homogeneous of degree two in W,X,Y,Z; no sphere reduction is required.",
       "psd_block_sizes":[b["order"] for b in blocks],
       "matrix_shape":[len(monomials),len(coordinates)],
       "monomial_rows":[list(p) for p in monomials],
       "gram_coordinates":coordinates,
       "entries":entries,
       "off_diagonal_convention":"A symmetric off-diagonal Gram coordinate contributes twice its feature product.",
       "coefficient_equation":"For a target with these normal-form coefficients b, A*vech(G)=b and every G block is positive semidefinite. All absent target coefficients are zero.",
    }
    unique_nonzero={tuple(sorted(f["polynomial"].items())) for f in original_features if f["polynomial"]}
    receipt={"original_named_features":len(original_features),"original_nonzero_features":sum(bool(f["polynomial"]) for f in original_features),
       "original_distinct_nonzero_normal_forms":len(unique_nonzero),
       "independent_fixed_functions":len(features),
       "polynomial_functions":sum("ijk" in f for f in features),
       "angular_functions":sum("n" in f for f in features),
       "psd_blocks":len(blocks),"total_order":sum(b["order"] for b in blocks),
       "block_size_counts":{str(n):sum(b["order"]==n for b in blocks) for n in sorted(set(b["order"] for b in blocks))},
       "maximum_block_order":max(b["order"] for b in blocks),
       "gram_coordinates":len(coordinates),"coefficient_rows":len(monomials),"nonzero_integer_entries":len(entries),
       "old_order669_recomputed":False,
       "all_products_evaluated_exactly":True,
       "not_a_Q222_pencil":True}
    for name,obj in [("WEIGHTED_FIXED_FUNCTIONS.json",data),("WEIGHTED_GRAM_PRODUCT_MAP.coo.json",product),("compilation_receipt.json",receipt)]:
        (HERE/name).write_text(json.dumps(obj,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
