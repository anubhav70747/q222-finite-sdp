"""Independent exact verifier for the order-33 separator.

Uses only Python's standard library and Fraction.  It does not use SymPy,
an SDP solver, the scout's objective implementation, or numerical signs.
"""
from fractions import Fraction as F
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("certificate", nargs="?", type=Path,
                    default=HERE / "dyadic_order33_separator.json")
PATH = parser.parse_args().certificate
cert = json.loads(PATH.read_text())

def mat(rows):
    return [[F(x) for x in row] for row in rows]

def tr(A):
    return sum(A[i][i] for i in range(len(A)))

def block(U, R, W):
    return [U[i]+R[i] for i in range(4)] + [R[i]+W[i] for i in range(4)]

def add(*As):
    return [[sum(A[i][j] for A in As) for j in range(len(As[0][0]))]
            for i in range(len(As[0]))]

def part(A, i, j):
    return [row[j:j+2] for row in A[i:i+2]]

def transpose(A):
    return list(map(list, zip(*A)))

def negate(A):
    return [[-x for x in row] for row in A]

hatted = [[mat(A) for A in row] for row in cert["hatted_moments"]]
original = [[mat(A) for A in row] for row in cert["original_moments"]]
scale = F(cert["bob_scale"])
diagonal = [[F(x) for x in row] for row in cert["state_congruences"]]
assert scale > 0 and all(x != 0 for row in diagonal for x in row)

for side in range(2):
    for k in range(3):
        H, O = hatted[side][k], original[side][k]
        assert H == transpose(H) and O == transpose(O)
        expected = [[scale**k*diagonal[side][i]*H[i][j]*diagonal[side][j]
                     for j in range(4)] for i in range(4)]
        assert O == expected

# Every PSD sign is proved by exact supplied rational LDL reconstruction.
min_pivots = []
for side in range(2):
    for name, A in [("M", block(*hatted[side])), ("R", hatted[side][1])]:
        proof = cert["hatted_LDL"][side][name]
        L = mat(proof["L"])
        p = [F(x) for x in proof["diagonal"]]
        n = len(A)
        assert all(pivot > 0 for pivot in p)
        if cert.get("rounding_denominator") == str(2**30):
            assert all(pivot > F(1, 2**27) for pivot in p)
        assert all(L[i][i] == 1 for i in range(n))
        assert all(L[i][j] == 0 for i in range(n) for j in range(i+1,n))
        reconstructed = [[sum(L[i][k]*p[k]*L[j][k] for k in range(n))
                          for j in range(n)] for i in range(n)]
        assert reconstructed == A
        min_pivots.append(("GH"[side]+name, min(p)))

# The twelve J equations, using all entries of each cross block.
for k in range(3):
    for i in range(2):
        for j in range(2,4):
            assert original[0][k][i][j] == original[1][k][i][j]

# Source coefficients computed by the independently derived block formulas.
S, V, C, D = [], [], [], []
for k in range(3):
    G, H = original[0][k], original[1][k]
    ag, bg, ah, bh = part(G,0,0), part(G,2,2), part(H,0,0), part(H,2,2)
    X = part(G,0,2)
    S.append(add(ag,ah))
    V.append(add(bg,bh))
    C.append(add(ah,bh,negate(X),negate(transpose(X))))
    D.append(add(ag,bg,X,transpose(X)))
y = [
    2*(S[0][1][1]+S[2][1][1]+V[0][0][0]+2*V[1][0][1]+V[2][1][1]
       +C[2][0][0]-2*C[1][0][1]-C[2][1][1]),
    4*(D[0][0][0]+2*D[1][0][1]+D[2][1][1]),
    4*(D[0][1][1]+D[2][1][1]),
    4*(S[2][0][0]-2*S[1][0][1]+S[0][1][1]),
    4*(V[0][0][0]+V[2][0][0]),
    4*(C[0][0][0]+C[2][0][0]),
    4*(V[0][0][0]+2*V[1][0][1]+V[2][1][1]),
    4*(C[2][0][0]-2*C[1][0][1]+C[0][1][1]),
    4*(S[0][1][1]+S[2][1][1])
]
assert y == [F(x) for x in cert["dual_source_coefficients"]]
assert sum(y[1:])/4-y[0] == sum(tr(row[0])+tr(row[2]) for row in original)
c = [F(x) for x in cert["source"]]
assert len(c) == 9 and all(x >= 0 for x in c)
assert cert["source_order"] == ["d", "lambda", "n", "m", "C", "A", "D", "alpha", "beta"]
assert c == [F(1), F(0), F(0), F(1166523000,31661820415017653),
             F(322803000,186199784281), F(218288946200,1359122513),
             F(90943434293,222507771414), F(248802430236,3300726103),
             F(205476973,51646655494)]
pairing = sum(a*b for a,b in zip(c,y))
assert pairing == F(cert["source_pairing"]) < 0

print("PASS: exact LDL reconstructions of four positive definite blocks.")
print("PASS: congruence conversion and all twelve selector overlap equations.")
print("PASS: independent nine-coordinate dual map and source signs.")
print("Exact negative pairing:", pairing)
print("Approximate pairing (display only):", float(pairing))

if cert.get("rounding_denominator") == str(2**30):
    # A simpler source in the strict positive orthant.  Coordinatewise
    # domination preserves validity because every added d=0 term is an
    # elementary nonnegative event square in the established canonical slack.
    cbar = [F(1), F(1,2**34), F(1,2**14), F(1,2**24), F(1,576),
            F(161), F(1,2), F(76), F(1,250)]
    assert all(a > 0 and a >= b for a,b in zip(cbar,c))
    ell = [2**42*a for a in y]
    assert all(a.denominator == 1 for a in ell)
    assert ell == [-23803618, 906393494300170, 2744379904,
                   1378590360912, 5935589516, 20425, 2600028,
                   91248, 278104928]
    cbar_pairing = sum(a*b for a,b in zip(ell,cbar))
    assert cbar_pairing == F(-5417798011632274183,9663676416000) < -560635
    print("PASS: all-positive source dominates the established valid source.")
    print("Exact all-positive-source pairing with the integer functional:", cbar_pairing)

    expected_floors = [
        [3376576,82695,48553,117468710,16,14,14,31],
        [2570,5564,1593,10646],
        [142545776,16,32,17,29078313,12,24,11],
        [201343592,8,15,8]
    ]
    for idx, (side,name) in enumerate((s,n) for s in range(2) for n in ("M","R")):
        scaled = [2**30*F(a) for a in cert["hatted_LDL"][side][name]["diagonal"]]
        floors = [a.numerator//a.denominator for a in scaled]
        equality_indices = [i+1 for i,a in enumerate(scaled) if a.denominator == 1]
        assert floors == expected_floors[idx]
        assert equality_indices == [1]
        print("LDL scaled-pivot floors", "GH"[side]+name, floors,
              "; lower-bound equality only at pivot 1")
