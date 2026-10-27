# C1 H1 global admission: exact scalar-polynomial gate and unresolved universal PSD

27 October 2026. **Not a proof of H1.** The preceding source-wide contact-jet theorem is preserved; the present work does not claim a universal degree-four skew, a K364 Gram, or a K364 separator. It makes the *one missing step of the particular 40-feature cross-parity mechanism* completely explicit in scalar source-coefficient polynomials and proves the endpoint algebra. Direct membership in K364 is weaker than this skew requirement.

## Physical source domain and full four-state parent

Take the native oriented C1 parameterization `k,f,g,a,t` with `0<k<=1`, `k,f,g,a,t>0`, `r=1-af^2>0`, `0<tg^2<1`, and the strict original inequalities `Hc,Qc,M,J>0`, plus the other genuine native positive denominator conditions of SD1/NC1. The six native retained event weights `m,C,A,D,alpha,beta` and multiplier `d` are strictly positive, with `lambda=n=0`. Use their *original rational coquotient functions*, not seven free weights for admission. The complete SG2 parent is

    F=4m(vW-X+uvY-uZ)^2+4C(1+v^2)(uW-Y)^2
     +4A(1+u^2)(1+v^2)Y^2
     +(4D+2d)(uW+uvX-Y-vZ)^2
     +(4alpha+2d)(1+u^2)(vY-Z)^2
     +(4beta+2d)(1+v^2)(X+uZ)^2
     -2d(1+u^2)(1+v^2)Z^2.

Set `h=1+v^2`, and order state pairs `w=(W,X)` and `y=(Y,Z)`. In `F=psi^T(P0+u P1+u^2 P2)psi`, the complete coefficient matrices are

    P0w = [[4m v^2, -4m v], [-4m v, 4m+(4beta+2d)h]]
    P0y = [[4(A+C)h+4D+2d h+4alpha v^2,4(D-alpha)v],
           [4(D-alpha)v,4D v^2+4alpha]]
    P2w = [[4C h+4D+2d,(4D+2d)v],[(4D+2d)v,(4D+2d)v^2]]
    P2y = [[4A h+(4alpha+2d+4m)v^2,-(4alpha+2d+4m)v],
           [-(4alpha+2d+4m)v,4alpha+4beta h+2d+4m]]
    E=P1[w,y]/2 = [[2m v^2-2C h-2D-d, -v(2D+d+2m)],
                  [-v(2D+d+2m), 2beta h-2D v^2+d+2m]].

For ANY source-only polynomial `T(v)` with 2x2 coefficients and degree at most 4, put `A0=h P0w`, `B0=h P0y`, `A2=h P2y`, `B2=h P2w`, `L=h E`, and the complete parity blocks

    HA(v) = [[A0,L+T],[(L+T)^T,A2]],
    HB(v) = [[B0,(L-T)^T],[L-T,B2]].

Their direct sum, after a constant permutation, represents `h F` **coefficientwise** for every Alice angle `u`, Bob angle `v` and all four real state coordinates, since the two skew contributions cancel on `(psi,u psi)`. There is no restriction to contact-normal subspaces.

## Exact positive denominators of both full parity Schur tests

The four physical source coefficients `m,beta,d,C,D,A,alpha` are positive on the strict domain. Direct expansion gives

    det P0w = 8 m(2 beta+d) v^2 h > 0,
    det P0y = 8 h [ 2(A+C)(D v^2+alpha) + 2 D alpha h
                     + d(D v^2+alpha)] > 0,
    det P2w = 8 C(2D+d) v^2 h > 0,
    det P2y = 8 h [2 A(alpha+beta h+m)+A d
                     +2 alpha beta v^2 + beta d v^2+2 beta m v^2] >0.

Thus `A0` and `B0` are positive definite for every **strictly positive** Bob angle, and every inverse used in the Schur analysis is an ordinary matrix inverse at the testing point. This inverse is not proposed as a Gram coefficient or mixed SOS denominator.

## Four source-polynomial inequalities equivalent to BOTH 2x2 Schur residuals

Put `CA=L+T`, `CB=L-T`. Define four scalar polynomials of Bob angle `v`, with coefficients algebraic in the **source only**:

    qA = det(A0)*tr(A2) - tr(CA^T adj(A0) CA),
    dA = det([[A0,CA],[CA^T,A2]]),
    qB = det(B0)*tr(B2) - tr(CB adj(B0) CB^T),
    dB = det([[B0,CB^T],[CB,B2]]).

Then, for every v>0, the **exact equivalence** is

    HA(v)>=0 and HB(v)>=0
        iff qA(v)>=0 and dA(v)>=0 and qB(v)>=0 and dB(v)>=0.

Proof: `RA=A2-CA^T A0^-1 CA`, `RB=B2-CB B0^-1 CB^T`. The two strict positive 2x2 leading blocks give Schur equivalence. Since `det(A0),det(B0)>0`, `tr(RA)=qA/det(A0)`, `det(RA)=dA/det(A0)`, with analogous B identities. A real symmetric 2x2 matrix is PSD iff its trace and determinant are nonnegative. Thus the four inequalities are individually necessary **and jointly sufficient** for the 40-feature skew construction. `qA,qB` have degree <=12, `dA,dB` degree <=16. These are polynomial inequalities, not merely positive *evaluated* Bell matrices.

For a half-line proof one may certify each scalar polynomial via ordinary one-variable Markov-Lukacs squares (`p=sigma0(v)+v sigma1(v)`), with degrees bounded by the displayed finite degrees. But no source-wide choice of the 20 skew coefficients and corresponding scalar SOS factors has been derived.

## Exact forced endpoint and moving-contact polynomial factors

If `HA(0),HB(0)` are PSD, the zero diagonal on `W0` in `HA(0)` and on `X1` in `HB(0)` forces their corresponding complete rows to vanish. This fixes the entire constant skew matrix:

    T(0) = diag(2C+2D+d, 2beta+d+2m).

For any T obeying this identity, both `dA(v)` and `dB(v)` are divisible by `v^2`: the unique zero row at v=0 has diagonal of order `v^2` and other entries of order `v`. If the same T matches the two actual native contact kernels at `v=v_*` and sets both first soft Bob derivatives to zero, each `dA,dB` has a double zero at `v_*`. Hence by exact monic polynomial division over the source field,

    dA(v)=v^2 (v-v_*)^2 pA(v),
    dB(v)=v^2 (v-v_*)^2 pB(v),

for source-polynomial `pA,pB` of degrees at most 12. This is an **algebraic consequence** of the two contact-jet and endpoint hypotheses, not a proof that `pA,pB>=0` between the proved neighbourhoods. If both traces `qA,qB` and the reduced determinants `pA,pB` are nonnegative on the whole positive half-line, then the two complete parity Grams are PSD globally and the independently verified degree-40 univariate matrix-SOS implication would prove **H1**. Those four scalar inequalities remain **unproved for any universal source-only T**.

## Exact statement that remains open

For every genuinely physical strict oriented native C1 source (`0<k<=1`), find a SINGLE source-only polynomial matrix `T_theta(v)` of degree <=4 such that `qA,dA,qB,dB` above are nonnegative for **all** `v>0` and the two homogeneous endpoints are PSD. Contact/endpoint choices may be made jointly, and all singular PSD certificates are permitted. The original source-wide contact-jet theorem does NOT imply these intermediate-angle polynomial inequalities.

If no such degree-four cross-parity skew exists, it disproves only this **sufficient 40-feature mechanism**, not H1. The alternative direct H1 assertion requires a 364x364 PSD Gram with source-only entries and the exact K364 coefficient map; no such universal Gram, and no K364 separating source, is obtained here. The fixed conditional outer SDP, existing qualitative existence proof and strict C2 open-end admission remain intact.

`check_scalar_gate.py` replays exactly the four endpoint-block determinant identities, the forced endpoint skew and the generic Schur/trace/determinant polynomial identities. It performs no parameter sweep, no SDP, and no test of an unproved global candidate.
