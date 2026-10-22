# Full original-source C1D allocation and a finite local certificate

28 September 2026. This is a complete result on the printed local domain, not a
whole-C1, whole-C2, or full-Q222 representation.

## 1. Original data, variables and scope

The unmodified inputs are `sources/XSMALL_Z_FULL_SOURCE.json`,
`sources/XSMALL_W_FULL_SOURCE.json`,
`sources/X_BOUNDED_V_SMALL_FULL_POLYNOMIALS.json`, and
`sources/ORIGINAL_C1_CONTACT.json`. The last two bind the source normalization to
the original physical C1 coefficient order `(d,e2,e3,e4,e5,alpha,beta)`.

Coefficient-side variables are `(a,b,q,r,E)`. Testing-side variables are
`(u,v,W,X,Y)`, with the independent tested state `(W,X,Y,1)`. In particular the
contact angle q is NOT identified with the tested angle v. Put h=a b. The source
contact substitution is

    u_contact=a b^2 q, v_contact=q, A_contact=1/a,
    Y_contact=r/a, e_contact=q^2 E.

Use 0<a<=1 on this original contact chart. The polynomial identity below also
extends to a=0 and negative a, but no claim about the contact parametrization is
made there. Coefficient nonnegativity is retained when identifying a point with
the physical master. A strictly positive original coefficient example inside
the chart is recorded in `data/SOURCE_ALLOCATION_VERIFICATION.json`.

## 2. Actual original-generator membership certificates

Write

    g=u W-b^2 v,  k=b u v,  s=b^2(v^2-q^2),
    I=(g,k,s),  n=(u,W,v,q,b),
    B=X+u-a^2 b^2 q^2 E v W,  N=Y-r u W.

The full source normal form has six arrays, R0,LX,LY,Rxy,HB,HN, and is

    Fhat=R0+LX B+LY N+Rxy B N
         +HB b^2 B^2/E + HN N^2/[r(1-r)].

This exhausts the full state-quadratic parent. No higher-order remainder is
omitted. The new original-generator identities are

    R0-8(g^2+k^2+s^2) in n I^2,
    LX in I^2+b n I,
    LY in I^2+n I.

For each target the release prints every polynomial multiplier against the
ORIGINAL generators of the stated ideal. The coefficient variables a,r,E remain
polynomial parameters. No nonconstant denominator is used by this computation.

| Target | Target terms | Parameter monomials | Original-multiplier terms | Max multiplier degree |
|---|---:|---:|---:|---:|
| R0-8(g^2+k^2+s^2) | 2045 | 107 | 13805 | 57 |
| LX | 568 | 71 | 6292 | 49 |
| LY | 884 | 70 | 6724 | 49 |

A labelled Buchberger producer constructs the identities. The independent
reader does NOT rely on completeness of its basis or a saved zero remainder:
it multiplies the printed original generators by their multipliers and compares
the result with each complete source target. All 26821 multiplier terms are
therefore verifiable without repeating the basis search.

## 3. Complete five-root Gram

Let e^2=E, rr^2=r, ss^2=1-r, with positive roots, and set

    z=(g,k,s,bB/e,N/(rr ss)).

The allocation produces the full rational polynomial matrix Q satisfying

    Fhat=z^T Q z.

The fifteen upper-triangle entries are printed in
`data/C1D_Z_COMPLETE_SOURCE_GRAM.json`. They contain 31028 terms in the explicit
13-variable ring `(u,W,v,q,b,a,r,E,X,Y,e,rr,ss)`, of maximum degree 58. This size
is an uncompressed exact output, not a minimality statement.

Top soft entries receive the scalar allocation and B or N times the corresponding
I^2 normal allocations. The b n I part of LX gives soft-to-hardB entries multiplied
by e/2; the n I part of LY gives soft-to-hardN entries multiplied by rr ss/2. The
mixed source array is divisible by b and gives

    Q_34=(Rxy/b)*e*rr*ss/2.

The two hard diagonals are the exact HB and HN source polynomials. Those hard
arrays are NOT replaced by their residues.

`code/verify_source_allocations.py` clears only the pure coefficient factors by
multiplying every root by e rr ss, reduces the three pure square-root relations,
and compares the result with an independently built original Pauli operator.
The resulting 5011-term polynomial identity passes exactly. This verification
checks the whole physical operator, not merely the three ideal memberships.

## 4. An explicit finite local box

Every nonconstant term of Q-8I contains at least one variable from

    (u,W,v,q,b,X,Y).

The sum of absolute rational coefficients over the upper triangle of Q-8I is

    T=1394008.

Choose

    epsilon=1/[2(1+T)]=1/2788018.

On the pure coefficient rectangle

    0<a<=1, 0<b,q<=epsilon, 0<r<1, 0<E<=1,
    e=sqrt(E), rr=sqrt(r), ss=sqrt(1-r),

and the independent testing rectangle

    |u|,|v|,|W|,|X|,|Y|<=epsilon, Z=1,

all other monomials have absolute value at most one. Each small-containing
monomial has absolute value at most epsilon. Each row sum of absolute corrections
is at most epsilon T<1/2. Therefore

    Q > (15/2) I_5.

This is a uniform finite inequality, not just a positive-residue assertion.
There is NO positive lower bound on E, r or 1-r. Their possibly arbitrarily weak
scales occur in the pure coefficient factors of z, while Q stays uniformly
positive. Literal E=0, r=0 and r=1 are outside this root chart and require their
own routing; they have not been added by a limiting SOS argument.

## 5. Positive pure coefficient restoration to the original support

Let D be half of the source's cleared e3 coefficient. The original bridge is
verified for all seven coefficients, and the exact factorization is

    D=r (1+b^2 r)(1+a^2 b^4 q^2)(b^2 q^2 r-1)
       (1-E a^2 b^2 q^2+a^2 b^4 q^2 r)
       (1-a^2 b^2 q^2+E a^2 b^2 q^4)
       (-1+b^2(1-r)-a^2 b^4 q^2(1-E q^2)).

It is positive on the displayed domain: the third and last factors are negative,
and all other factors are positive. In particular the last factor is below
-1+b^2<0; the two preceding factors are bounded below by 1-a^2 b^2 q^2>0.

The independently verified parent relation is

    Fhat = D/[r(1-r)] * psi^T P_original(u,v) psi.

Hence

    psi^T M_original psi
      = [r(1-r)/D] * z^T Q z / [(1+u^2)(1+v^2)].

The first normalization is pure coefficient-side; the second is pure testing-side.
For the normalized physical state, also divide by 1+W^2+X^2+Y^2; its square
root is purely testing-side and is included in the normalized local testing
algebra export. Neither restoration is a mixed inverse. No dependence of final testing functions on coefficient
values is introduced by this restoration.

## 6. What this does and does not close

The release extracts the actual full-parent coefficient allocation, a uniform
finite domain, all pure-side square factors, and a fixed local testing algebra.
The identity gives a full local coefficient/testing tensor certificate for every
point in the printed rectangle. It is not a certificate for all independent
measurements or states at a fixed coefficient point. The finite union of testing
charts and the finite support-family selection required for global compilation
remain to be constructed and proved exhaustive.
