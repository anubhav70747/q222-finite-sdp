#!/usr/bin/env python3
"""Generate an UNSOLVED, exact first-order input for a three-block SOC lift.

No quantifier elimination is performed. Standard-library exact arithmetic.
Fresh checks verify the Pauli atom map, degree/count data and deterministic
behavior affine rank. They do not verify feasibility of the candidate lift.
"""
import hashlib
import itertools
import json
from pathlib import Path
from exact_poly import P, expectation, kron, mat_add, mat_scale, rank_rational

HERE = Path(__file__).resolve().parent
COORDINATES = ["a0", "a1", "b0", "b1", "E00", "E01", "E10", "E11"]
checks = []


def check(name, statement):
    if not statement:
        raise AssertionError(name)
    checks.append({"name": name, "status": "PASS"})


def write(name, obj):
    (HERE/name).write_text(json.dumps(obj, indent=2, sort_keys=True)+"\n")


def atom(names):
    w, x, y, z, ca, sa, cb, sb = map(P.var, names)
    ta, xa = w*w+x*x-y*y-z*z, 2*(w*y+x*z)
    tb, xb = w*w-x*x+y*y-z*z, 2*(w*x+y*z)
    czz, czx = w*w-x*x-y*y+z*z, 2*(w*x-y*z)
    cxz, cxx = 2*(w*y-x*z), 2*(w*z+x*y)
    phi = [ta, ca*ta+sa*xa, tb, cb*tb+sb*xb,
           czz, cb*czz+sb*czx, ca*czz+sa*cxz,
           ca*cb*czz+ca*sb*czx+sa*cb*cxz+sa*sb*cxx]
    constraints = [w*w+x*x+y*y+z*z-1, ca*ca+sa*sa-1, cb*cb+sb*sb-1]
    return phi, constraints


local_names = ["w", "x", "y", "z", "ca", "sa", "cb", "sb"]
phi, norms = atom(local_names)
w, x, y, z, ca, sa, cb, sb = map(P.var, local_names)
Z, X, I = [[1, 0], [0, -1]], [[0, 1], [1, 0]], [[1, 0], [0, 1]]
A1, B1 = mat_add(mat_scale(ca, Z), mat_scale(sa, X)), mat_add(mat_scale(cb, Z), mat_scale(sb, X))
operators = [kron(Z, I), kron(A1, I), kron(I, Z), kron(I, B1),
             kron(Z, Z), kron(Z, B1), kron(A1, Z), kron(A1, B1)]
for name, op, target in zip(COORDINATES, operators, phi):
    check("pauli_atom_"+name, expectation([w, x, y, z], op) == target)
check("atom_degree_vector", [p.degree for p in phi] == [2, 3, 2, 3, 2, 3, 3, 4])

qnames = ["q222"+name for name in COORDINATES]
qvars = list(map(P.var, qnames))
weights = ["q222lambda"+str(i) for i in range(1, 10)]
atom_names = [["q222"+name+str(i) for name in local_names] for i in range(1, 10)]
physical_vars = weights + [name for group in atom_names for name in group]
atoms = [atom(names) for names in atom_names]
equalities = [sum((P.var(name) for name in weights), P())-1]
equalities += [p for _, constraints in atoms for p in constraints]
for index, q in enumerate(qvars):
    equalities.append(q-sum((P.var(weight)*group[0][index]
                            for weight, group in zip(weights, atoms)), P()))
inequalities = list(map(P.var, weights))
check("physical_quantified_variable_count", len(physical_vars) == len(set(physical_vars)) == 81)
check("physical_equality_count", len(equalities) == 36)
check("physical_inequality_count", len(inequalities) == 9)
check("physical_maximum_degree", max(p.degree for p in equalities+inequalities) == 5)
check("physical_output_degree_vector", [p.degree for p in equalities[-8:]] == [3, 4, 3, 4, 3, 4, 4, 5])

# A rational affine-rank calculation, separate from the quantum atom map.
deterministic = []
for av0, av1, bv0, bv1 in itertools.product([-1, 1], repeat=4):
    deterministic.append([av0, av1, bv0, bv1, av0*bv0, av0*bv1, av1*bv0, av1*bv1])
affine_rank = rank_rational([[1]+point for point in deterministic])-1
check("full_quantum_affine_dimension", affine_rank == 8)

M = 3
D = 3*M
aux_names = ["q222aux"+str(i) for i in range(1, D+1)]
affine_basis = [P(1)] + qvars + list(map(P.var, aux_names))
coefficient_count = 3*M*(9+D)
coefficient_names = ["q222c"+str(i).zfill(3) for i in range(1, coefficient_count+1)]
blocks, coefficient_map = [], []
index = 0
for block in range(M):
    entries = []
    for entry in ["aa", "bb", "cc"]:
        names = coefficient_names[index:index+len(affine_basis)]
        value = sum((P.var(name)*term for name, term in zip(names, affine_basis)), P())
        entries.append(value)
        coefficient_map.append({"block": block+1, "entry": entry, "coefficient_symbols": names})
        index += len(affine_basis)
    blocks.append(entries)
check("soc_coefficient_count", index == coefficient_count == 162)
check("soc_auxiliary_count", len(aux_names) == 9)
check("soc_affine_basis_count", len(affine_basis) == 18)
check("all_variable_names_disjoint", len(set(physical_vars+qnames+aux_names+coefficient_names))
      == len(physical_vars+qnames+aux_names+coefficient_names))
psd = [constraint for aa, bb, cc in blocks for constraint in [aa, cc, aa*cc-bb*bb]]
check("soc_psd_scalar_constraints", len(psd) == 9)

def wl_list(items):
    return "{"+", ".join(items)+"}"


physical_body = "And[\n  "+",\n  ".join([p.wolfram()+" == 0" for p in equalities]
                                              +[p.wolfram()+" >= 0" for p in inequalities])+"\n]"
lift_body = "And[\n  "+",\n  ".join(p.wolfram()+" >= 0" for p in psd)+"\n]"
wl = """(* Q222_SOC_M3_UNSOLVED_INPUT.wl
   Regenerated 2026-10-22. Exact rational input; NO solver has been run.
   Three real symmetric 2x2 PSD blocks, 162 coefficient unknowns, 9 auxiliaries.
   Physical membership uses 9 real two-qubit atoms and 81 quantified variables.
   Loading this file only constructs formulas. It does not execute Resolve,
   Reduce, FindInstance or any optimization. Global extraction remains OPEN.
*)

"""
physical_formula = "Exists["+wl_list(physical_vars)+",\n"+physical_body+"]"
lift_formula = "Exists["+wl_list(aux_names)+",\n"+lift_body+"]"
wl += "q222PhysicalFormula = "+physical_formula+";\n\n"
wl += "q222LiftFormula = "+lift_formula+";\n\n"
# Expand quantifier bodies literally: lexical substitution through nested
# scoping constructs could otherwise rename a binder and change capture.
wl += "q222CandidateSentence = Exists["+wl_list(coefficient_names)+",\n"
wl += "  ForAll["+wl_list(qnames)+",\n    Equivalent[\n"
wl += physical_formula+",\n"+lift_formula+"\n]]];\n\n"
wl += "(* Deliberately inactive: an exact decision call could be\n"
wl += "   Resolve[q222CandidateSentence, Reals]\n"
wl += "   Even a True decision would still require extracting exact coefficients\n"
wl += "   and independently verifying the resulting pencil before delivery.\n*)\n"
(HERE/"Q222_SOC_M3_UNSOLVED_INPUT.wl").write_text(wl)
spec = {
    "status": "UNSOLVED", "quantifier_elimination_run": False, "coefficient_tuple_extracted": False,
    "global_E1_E2_E3": "OPEN", "purpose": "Exact first-order specification of a candidate three-block SOC representation of full Q222",
    "existence_premise": "Not freshly verified. M=3 feasibility is neither assumed nor established.",
    "coordinate_order": COORDINATES, "coordinate_symbols": qnames,
    "physical_membership": {
        "quantifier": "exists", "variables": physical_vars, "variable_count": len(physical_vars),
        "real_two_qubit_atom_count": 9, "weights": weights, "atom_variable_groups": atom_names,
        "local_atom_order": local_names, "atom_map": [p.json_terms() for p in phi],
        "equalities_zero": [p.json_terms() for p in equalities],
        "inequalities_nonnegative": [p.json_terms() for p in inequalities],
        "equality_count": len(equalities), "inequality_count": len(inequalities), "maximum_total_degree": 5},
    "candidate_lift": {
        "block_orders": [2]*M, "matrix_order_total": 2*M, "matrix_space_dimension": D,
        "auxiliary_symbols": aux_names, "auxiliary_count": len(aux_names),
        "coefficient_symbols": coefficient_names, "coefficient_count": coefficient_count,
        "affine_basis": ["1"]+qnames+aux_names, "coefficient_map": coefficient_map,
        "psd_block_convention": "[[aa,bb],[bb,cc]] is PSD iff aa>=0,cc>=0,aa*cc-bb^2>=0",
        "block_entry_polynomials": [[p.json_terms() for p in block] for block in blocks]},
    "sentence": "exists coefficients: forall q: (exists physical variables: physical constraints) iff (exists auxiliaries: all blocks PSD)",
    "scope_of_fresh_checks": "Algebraic encoding, atom expectations, exact counts/degrees and affine rank only; not feasibility or quantifier elimination",
}
write("Q222_SOC_M3_UNSOLVED_INPUT.json", spec)
write("CHECK_Q222_QE_ENCODING.json", {
    "status": "PASS_FOR_ENCODING_ONLY", "quantifier_elimination_run": False,
    "candidate_feasibility": "UNSOLVED", "exact_coefficients": "NOT_EXTRACTED",
    "check_count": len(checks), "checks": checks,
    "physical_formula_counts": {"variables": 81, "equalities": 36, "inequalities": 9, "max_degree": 5},
    "candidate_counts": {"blocks": 3, "coefficient_unknowns": 162, "auxiliary_variables": 9},
    "full_affine_dimension": affine_rank,
})
names = ["exact_poly.py", "generate_q222_qe_problem.py", "Q222_SOC_M3_UNSOLVED_INPUT.wl",
         "Q222_SOC_M3_UNSOLVED_INPUT.json", "CHECK_Q222_QE_ENCODING.json"]
write("QE_ENCODING_MANIFEST.json", {"schema": "sha256-v1", "status": "UNSOLVED_INPUT_ONLY", "files": {
    name: {"sha256": hashlib.sha256((HERE/name).read_bytes()).hexdigest(), "bytes": (HERE/name).stat().st_size}
    for name in names}})
print(json.dumps({"status": "PASS_FOR_ENCODING_ONLY", "check_count": len(checks),
                  "physical_variables": 81, "equalities": 36, "inequalities": 9,
                  "maximum_degree": 5, "candidate_coefficients": 162, "candidate_auxiliaries": 9,
                  "qe_status": "UNSOLVED"}, sort_keys=True))
