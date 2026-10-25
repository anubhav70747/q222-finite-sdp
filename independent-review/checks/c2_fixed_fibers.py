"""Exact fixed-(a,k) contact-fiber checks; no solver or private imports.

Requires SymPy. This verifies identities supporting the accompanying proof,
not the full three-parameter convex hull or a global SDP lift.
"""
import json
from pathlib import Path
import sympy as S


def main():
    checks = {}

    def zero(expr, name):
        assert S.cancel(expr) == 0, name
        checks[name] = True

    a, r, s = S.symbols("a r s", positive=True)
    e, p = r+s, r*s
    g = a+e-1
    L = a*(e-p)+p
    Q = a*g+L
    Delta = a**3*p+a*a*(e*e-p)+a*p*(2*e-p)+p*p
    B0 = a**3*p+a*a*(e*e+p-e)+a*((3*e-2)*p-p*p)+p*p
    native = [S.Integer(1), -a*B0/(2*(a+1)*Delta),
              (a+s)*(1-r)*(a*s+r)*L/(2*(a+1)*g*Delta),
              r*r*(a+s)*Q/(2*(1-r)*Delta),
              a*(a+r)*(a+s)*L/(2*(a+1)*Delta),
              (a+r)*(1-s)*(a*r+s)*L/(2*(a+1)*g*Delta),
              s*s*(a+r)*Q/(2*(1-s)*Delta)]
    xn, yn = (a+r)/(1-r), (a+s)/(1-s)
    kn = xn*yn
    d, n, m, C, A, D, al = native
    zero(a*(d+2*n)-2*A, "native_angle")
    zero(C+d/2-(kn-1)*D, "native_C_ratio")
    zero(al+d/2-(kn-1)*m, "native_alpha_ratio")
    zero((kn-a*a)*(A+n)-a*a*(kn-1)*(m+D),
         "native_fourth_linear_relation")
    Xn = kn/(kn-a)*((a-1)*A+a*(kn-1)*D)
    Yn = kn/(kn-a)*((a-1)*A+a*(kn-1)*m)
    zero(Xn-A*xn, "native_X_equals_Ax")
    zero(Yn-A*yn, "native_Y_equals_Ay")

    x, y, k, w = S.symbols("x y k w", positive=True)
    inverse = {r: (x-a)/(x+1), s: (y-a)/(y+1)}
    J = (x*y-a)*(x+y)+2*(1-a)*x*y
    Phi = x*y*(x*y-a*a)-a*J
    zero(g.subs(inverse)-(a+1)*(x*y-1)/((x+1)*(y+1)),
         "g_coordinate_identity")
    zero(B0.subs(inverse)-(a+1)**3*Phi/((x+1)**2*(y+1)**2),
         "B0_coordinate_identity")
    wminus = k*(k+a*a-2*a)/(a*(k-a))
    wplus = a+k/a
    zero((k-a*a)/(k-a)-(wplus-wminus), "nonempty_interval_gap")
    zero(k+a*a-2*a-((a-1)**2+k-1), "wminus_positive_numerator")

    # Four equalities uniquely reconstruct a source from A,X,Y.
    AA, XX, YY = S.symbols("A X Y", real=True)
    DD = ((k-a)*XX/k-(a-1)*AA)/(a*(k-1))
    mm = ((k-a)*YY/k-(a-1)*AA)/(a*(k-1))
    nn = a*a*(k-1)*(mm+DD)/(k-a*a)-AA
    dd = 2*AA/a-2*nn
    CC, ll = (k-1)*DD-dd/2, (k-1)*mm-dd/2
    reconstructed = [dd, nn, mm, CC, AA, DD, ll]
    cw = a*(k-a)/(k*(k-a*a))
    zero(nn-cw*(XX+YY-AA*wminus), "n_lower_w_bound")
    zero(a*AA-k*nn-k*cw*(AA*wplus-XX-YY), "upper_w_bound")
    zero(dd-2*(k-a*a)*AA/(a*k)-2*(a*AA-k*nn)/k,
         "strict_positive_d_bound")
    for i, expr in enumerate(reconstructed):
        zero(expr.subs({XX: xn, YY: yn, AA: 1, k: kn})
             -native[i]/A, f"all_native_reconstruction_{i}")

    rank_one = {AA: 1, XX: x, YY: k/x}
    zero(CC.subs(rank_one)-(k-a)*(x-a)**2/(a*x*(k-a*a)),
         "rank_one_C_positive_formula")
    zero(ll.subs(rank_one)-(k-a)*(k/x-a)**2/(a*(k/x)*(k-a*a)),
         "rank_one_alpha_positive_formula")
    zero(DD.subs({AA: 1, XX: a})-(k-a*a)/(k*a*(k-1)),
         "D_at_smallest_x_positive")
    # At fixed w, the PSD determinant is the interval condition for X/A.
    lo, hi, v = S.symbols("lo hi v", real=True)
    zero((v-lo)*(hi-v)-(v*(lo+hi-v)-lo*hi),
         "two_generator_interval_identity")
    zero((w-2*S.sqrt(k))*(w+2*S.sqrt(k))-(w*w-4*k),
         "real_endpoint_discriminant")
    zero((a-x)*(a-y)-(a*a-a*(x+y)+x*y),
         "endpoint_lower_bound_identity")

    # Quartic-to-PSD conversion uses only the four linear equalities.
    BB, UU, VV, MM = AA+nn, CC+dd/2, ll+dd/2, mm+DD
    SS = UU+VV
    R = ((k-a*a)*AA**2+a*(1-a)*k*AA*MM
         -a*a*k*(k-1)*mm*DD)
    quartic = SS**2*(BB+dd/2)*(MM+dd/2)-(a+1)**2*UU*VV*MM*(BB+SS)
    zero(quartic-(a+1)**2*(k-1)**2*MM**2*R/(a*a*(k-a*a)),
         "quartic_factor_on_linear_fiber")
    zero(XX*YY-k*AA**2+k*(k-1)*R/(k-a)**2,
         "PSD_determinant_vs_quartic_R")

    # Two exact examples printed in the received answer.
    wrong = list(map(S.Rational, [1, 38, "77/4", "33/4", 77, "7/8", 192]))
    base = list(map(S.Rational, [1, "43/299", "77/1495", "891/598",
                                "385/299", "119/598", "9/598"]))
    def fiber_values(c, av, kv):
        di, ni, mi, ci, ai, Di, ali = c
        Xi = kv/(kv-av)*((av-1)*ai+av*(kv-1)*Di)
        Yi = kv/(kv-av)*((av-1)*ai+av*(kv-1)*mi)
        return [av*(di+2*ni)-2*ai, ci+di/2-(kv-1)*Di,
                ali+di/2-(kv-1)*mi,
                (kv-av**2)*(ai+ni)-av**2*(kv-1)*(mi+Di)], Xi, Yi
    relations, Xw, Yw = fiber_values(wrong, S.Integer(2), S.Integer(11))
    assert relations == [0]*4 and Xw*Yw == 11*wrong[4]**2
    assert 2*wrong[4]-11*wrong[1] == -264
    checks["wrong_branch_example"] = True
    swapped = [base[i] for i in (0, 1, 5, 6, 4, 2, 3)]
    midpoint = [(v+q)/2 for v, q in zip(base, swapped)]
    di, ni, mi, ci, ai, Di, ali = midpoint
    bi, ui, vi, Mi = ai+ni, ci+di/2, ali+di/2, mi+Di
    Si = ui+vi
    residual = S.factor(Si**2*(bi+di/2)*(Mi+di/2)-9*ui*vi*Mi*(bi+Si))
    assert residual == -S.Rational(1543042787517, 319701552040)
    relations, Xm, Ym = fiber_values(midpoint, S.Integer(2), S.Integer(11))
    assert relations == [0]*4 and Xm*Ym > 11*ai**2
    checks["midpoint_is_strict_PSD_and_not_quartic_zero"] = True

    out = {
        "checks": checks,
        "number_of_exact_checks": len(checks),
        "fixed_parameter_domain": "a>0; k>max(1,a^2)",
        "original_generator_domain": "x,y>a; xy=k; x+y>=wminus (hence x+y<wplus)",
        "closed_generator_domain": "x,y>=a; xy=k; x+y>=wminus",
        "wminus": str(wminus), "wplus": str(wplus),
        "n_over_A": "a*(k-a)*(w-wminus)/(k*(k-a^2))",
        "d_over_A_lower_bound": "2*(k-a^2)/(a*k)>0",
        "C_over_A_on_generator": "(k-a)*(x-a)^2/(a*x*(k-a^2))",
        "alpha_over_A_on_generator": "(k-a)*(y-a)^2/(a*y*(k-a^2))",
        "generator_count_for_each_nonzero_point": "at most two, permitting endpoint limits",
        "quartic_midpoint_residual": str(residual),
        "scope": "Exact fixed-fiber closed conic hull only; no finite lift of the continuum of fibers."
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
