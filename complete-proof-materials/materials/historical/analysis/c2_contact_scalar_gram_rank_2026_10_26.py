"""Three specified native contacts; exact contact-scalar Gram ranks.

No solver, sweep, or general five-parameter expansion. Requires SymPy.
The general rank-three proof is in the matching text note.
"""
from pathlib import Path
import json
import sympy as S


def main():
    u, v, s, t = S.symbols("u v s t", real=True)
    cases = [
        ("nonME", [1, "7/66", "5/33", "5/22", "10/33", "16/33", "9/22"],
         1, 1, ["5/4", "-1/4", "1/2", 1]),
        ("ME_equal_angles", [1, "1/16", "3/16", "3/4", "9/16", "7/16", "1/28"],
         "1/2", "1/2", [1, "-1/4", "1/4", 1]),
        ("ME_unequal_angles", [1, "43/299", "77/1495", "891/598", "385/299", "119/598", "9/598"],
         "1/4", "1/2", [1, "-1/8", "1/8", 1]),
    ]
    output = []
    for label, cv, uc, vc, pv in cases:
        d, n, m, C, A, D, al = map(S.Rational, cv)
        uc, vc = S.Rational(uc), S.Rational(vc)
        W, X, Y, Z = map(S.Rational, pv)
        delta = W*Z-X*Y
        Du, Dv = 1+u*u, 1+v*v
        F = (4*n*Du*Dv*X*X+4*m*(v*W-X+u*v*Y-u*Z)**2
             +4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y*Y
             +(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2
             +(4*al+2*d)*Du*(v*Y-Z)**2
             +2*d*Dv*(X+u*Z)**2-2*d*Du*Dv*Z*Z)
        poly = S.Poly(S.expand(F.subs({u: s+uc, v: t+vc})/(4*delta**2)), s, t)
        co = lambda i, j: poly.coeff_monomial(s**i*t**j)
        assert co(0, 0) == co(1, 0) == co(0, 1) == 0
        G = S.Matrix([[co(2, 0), co(1, 1)/2, co(2, 1)/2],
                      [co(1, 1)/2, co(0, 2), co(1, 2)/2],
                      [co(2, 1)/2, co(1, 2)/2, co(2, 2)]])
        q = S.Matrix([s, t, s*t])
        assert S.expand(poly.as_expr()-(q.T*G*q)[0]) == 0
        leading = [S.factor(G[:i, :i].det()) for i in (1, 2, 3)]
        assert all(value > 0 for value in leading)
        output.append({
            "case": label, "active_order": ["d", "n", "m", "C", "A", "D", "alpha"],
            "source": list(map(str, cv)), "u_contact": str(uc), "v_contact": str(vc),
            "state": list(map(str, pv)), "det_Psi": str(delta),
            "G_gamma": [[str(x) for x in G.row(i)] for i in range(3)],
            "positive_leading_minors": list(map(str, leading)), "rank": 3
        })
    result = {
        "examples": output,
        "all_three_exact_Gram_identities_and_strict_positive_minors_pass": True,
        "general_rank_result": "rank G_gamma=3 for every finite native C2 contact with strict gaps, including n=0; see the analytical proof in the matching .txt",
        "scope": "No uniform eigenvalue bound, complete state-coupling certificate, or multiplier removal."
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
