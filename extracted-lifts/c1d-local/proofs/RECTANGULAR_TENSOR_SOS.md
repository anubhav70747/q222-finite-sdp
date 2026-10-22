# From a bounded source Gram to explicit pure-side squares

## 1. The matrix decomposition

For the five-root source Gram Q, use the exact identity

    Q = sum_i (Q_ii-4) e_i e_i^T
        + sum_(i<j) [(1+Q_ij)/2] (e_i+e_j)(e_i+e_j)^T
        + sum_(i<j) [(1-Q_ij)/2] (e_i-e_j)(e_i-e_j)^T.

There are 25 scalar channels. Their positive reference values are four for a
diagonal channel and one-half for an off-diagonal channel. The full finite-box
bound from EXPLICIT_C1D_SOURCE_ALLOCATION.md makes every scalar positive, but
ordinary positivity alone is not used to identify a mixed element as a square.

## 2. Explicit tensor splitting with no mixed radical

In a channel write

    q= rho + sum_m beta_m C_m tensor T_m,
    S=sum_m |beta_m|.

C_m and T_m are the coefficient and testing monomials, respectively. Every term
has at least one small pure factor. If C_m contains b or q, use coefficient-small
orientation and put t=sign(beta_m) C_m, l=T_m. Otherwise T_m contains one of the
small testing variables and the factors are interchanged. In either case |t|<=epsilon,
|l|<=1, and the exact bound gives rho/S>epsilon.

For each term use

    rho/S + t tensor l
      = [(rho/S+t) tensor (1+l)
         +(rho/S-t) tensor (1-l)]/2.

If the testing side is small, the two displayed tensor factors are interchanged.
All four radicands are nonnegative IN THEIR OWN real closed factors. Multiplying
by |beta_m| and summing gives the original scalar channel. Taking the separate
square roots and multiplying by the corresponding root combination e_i^T z or
(e_i+/-e_j)^T z produces a tensor-SOS identity.

Each of the five z_i is explicitly a sum of at most two simple tensors. Thus
every resulting root has tensor separation rank at most four on this local
coefficient/testing domain. This local rank bound is not a global sxdeg bound.

## 3. Fully printed local library

`C1D_Z_EXPLICIT_SEPARATED_SQUARE_LIBRARY.json` prints the finite index set of 41464
term pairs, all rational weights, both monomial exponent lists, the small-side
assignment, scalar budget, root direction, and testing-radical references. Its
formula therefore specifies 82928 individual tensor squares, without an implicit
choice of a sufficiently large or small constant.

The independently reconstructed local testing algebra has:

* 1212 distinct testing radicals;
* 3844 fixed testing spanning functions;
* 1223 algebra variables and 1218 polynomial equations;
* maximum relation degree 23.

The feature functions are

    kappa * zeta_j * chi_a,

where chi_a is one of

    1, uW, v, uv, v^2, X+u, vW, Y,

and only the 3844 printed (j,a) pairs are retained. This is a spanning list, not
an asserted linearly independent basis. The normalizer obeys

    kappa^2(1+u^2)(1+v^2)(1+W^2+X^2+Y^2)=1.

This variant gives the slack of the normalized physical state. Omitting the last
factor gives the unnormalized state quadratic form. Both normalizers are pure
testing-side and have denominators uniformly bounded away from zero.

For each testing coordinate t, a box variable satisfies

    box_t^2+t^2=epsilon^2.

Each zeta_j^2 equals the printed nonnegative testing radicand. These equations
cover the entire closed local testing box, and all additional real sign choices
leave the square identity unchanged. The variety is compact: the box coordinates,
all radicands and the strictly positive normalizing denominator are bounded.

Some radicals depend on normalized state components as well as angles. This is
legal in the algebraic/testing-space route and is not restricted to the earlier
angle-only, state-linear matrix-SOS ansatz.

## 4. Independent verification

`verify_local_square_library.py` independently reconstructs the full five-by-five
source Gram from every printed pair, checks all rational scalar budgets, verifies
ownership and small-factor guards, checks the separated five-root expansion,
and exports the complete local testing algebra and spanning list.

The resulting identity is for the ORIGINAL unmultiplied physical slack after
one pure coefficient and one pure testing restoration. It is not an identity
for sigma3(M) times the slack, and no such mixed multiplier is cancelled.

## 5. Remaining selection problem

This library is common for the displayed LOCAL rectangle only. It is not yet a
common global library for all physical supports and all independent tests. The
other source charts need equivalent explicit data, followed by finite test-image
and support-admission selections with proved coverage. No positivity test on this
one rectangle can replace those operations.
