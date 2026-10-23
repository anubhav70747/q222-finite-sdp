"""Read-only exact audit of the already generated coefficient-bound receipts."""
from pathlib import Path
from fractions import Fraction as F
import json
base=Path(__file__).parent
data=base/"data"
u=json.loads((data/"lII_uniform_bound.json").read_text())
g=json.loads((data/"LII_FIXED_GATE_BUDGETS.json").read_text())
source=F(u["source_sum"]); testing=F(u["testing_sum"])
assert source>=0 and testing>=0
assert source<F(1,2**40)
assert all(F(r["source_bound"])>=0 and F(r["testing_bound"])>=0
           for r in u["records"])
expected_keys=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),
 (2,1),(2,2),(3,0),(3,2),(4,0),(4,2),(5,0)]
expected_exps=[40,50,33,4,55,41,18,62,48,28,65,42,72,51]
assert [tuple(r["exponents"]) for r in g["budgets"]]==expected_keys
for row,power in zip(g["budgets"],expected_exps):
    assert F(row["budget"])==F(1,2**power)
    assert 0<=F(row["raw_bound"])<=F(row["budget"])
total=sum(F(row["budget"]) for row in g["budgets"])
assert total==F(295165937727343559809,4722366482869645213696)
assert total<F(1,8)
out={
 "matrix_entry_count":len(u["records"]),
 "source_sum_exact":str(source),"source_sum_lt_2_minus_40":True,
 "testing_sum_exact":str(testing),
 "all_entry_bounds_nonnegative":True,
 "all_14_budgets_match_recovered_SC7":True,
 "total_budget_exact":str(total),"total_budget_lt_1_over_8":True,
 "reserve_lower_bound":str(F(1,2)-total),
 "reserve_strictly_above_3_over_8":True,
}
if (data/"LII_NORMALIZED_SOURCE_WITNESS.json").exists():
    w=json.loads((data/"LII_NORMALIZED_SOURCE_WITNESS.json").read_text())
    rows=w["upper_entries_by_test_monomial"]
    assert len(rows)==51
    assert sorted(map(tuple,w["testing_support"]))==expected_keys
    assert all(0<=i<=j<6 and tuple(e) in expected_keys
               for i,j,e,c in rows)
    out["literal_export_records"]=51
    out["literal_export_support"]=14
(base/"strong_bound_receipt_audit.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
