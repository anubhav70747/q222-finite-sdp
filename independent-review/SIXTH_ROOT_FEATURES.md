# Which fractional functions are needed?

The Motzkin limit explains why fixed NPA word lists fail. It does not
determine which fractional functions will certify the complete parent.
Two exact identities separate these questions.

## The scalar limit already admits half-powers

For positive $u,v$, the Motzkin polynomial satisfies

$$
\begin{aligned}
M(u,v)&=1+u^4v^2+u^2v^4-3u^2v^2\\
&=(u^2v-uv^2)^2
+2\bigl((uv)^{3/2}-(uv)^{1/2}\bigr)^2+(uv-1)^2.
\end{aligned}
\tag{1}
$$

A nonnegative constant adds one constant square. Thus the scalar limiting
polynomial, on the positive quadrant, cannot require thirds or sixths.
The fractional factors in (1) are not polynomial functions across the
commuting axes, so (1) does not supply the polynomial SOS excluded by the
NPA argument.

Nor does (1) place a complete Bell parent in the half-power cone. The
original testing span has angular exponents at most one and roots linear
in the state. The scalar roots in (1) reach angular exponent two.
Equation (1) describes the limiting polynomial obtained after the state
substitution and rescaling; it supplies no inverse construction of
admissible roots for the full parent. The [exact C2 half-power exclusion](C2_HALF_POWER_OBSTRUCTION.md)
therefore concerns an additional obstruction involving the full state matrix.

## The sixth-root grid has four weighted polynomial blocks

Put $p=u^{1/3}$ and $q=v^{1/3}$. The 196 functions in $Z_6$ form the four
groups

| Functions, for each state coordinate $\psi_\ell$ | Index range | Number |
|---|---|---:|
| $p^i q^j\psi_\ell$ | $0\leq i,j\leq3$ | 64 |
| $\sqrt p\,p^i q^j\psi_\ell$ | $0\leq i\leq2,\ 0\leq j\leq3$ | 48 |
| $\sqrt q\,p^i q^j\psi_\ell$ | $0\leq i\leq3,\ 0\leq j\leq2$ | 48 |
| $\sqrt{pq}\,p^i q^j\psi_\ell$ | $0\leq i,j\leq2$ | 36 |

After $u=r^6,v=s^6$, sign averaging separates these four parities.
Consequently the full Gram condition is exactly a weighted polynomial
SOS in the cubic coordinates, with weights $1,p,q,pq$ and the displayed
degree bounds. The three functions in one odd row are
$\sqrt p,p\sqrt p,p^2\sqrt p$, or $u^{1/6},u^{1/2},u^{5/6}$.
Retaining the middle function alone can lose certificates.

## A smaller list can lose a biquadratic certificate

Consider

$$P(u,v)=1+u^2v+uv^2-3uv\geq0\qquad(u,v\geq0).$$

This example is not asserted to be a native Bell source. It proves that
the sixth-root grid cannot automatically be replaced by the smaller
per-angle exponent set $S=\{0,1/3,1/2,2/3,1\}$ for arbitrary
state-linear biquadratic forms.

First, $PZ^2$ has a certificate in $Z_6$. Set
$A=r^4s^2$, $B=r^2s^4$, $C=1$, where $u=r^6,v=s^6$. The identity

$$
A^3+B^3+C^3-3ABC
=\frac12\sum_{h\in\{A,B,C\}}\sum_{i<j\in\{A,B,C\}}
\bigl[\sqrt h\,(i-j)\bigr]^2
\tag{2}
$$

gives nine roots. All their monomials have separate exponents between
zero and six in $r,s$. Multiplying every root by the state coordinate
$Z$ therefore gives the required state-linear certificate.

Now suppose an SOS used only $S\times S$. In $r,s$ its root exponents
would belong to $J^2$, where $J=\{0,2,3,4,6\}$. Every root exponent
lies in half the Newton polygon of
$1+r^{12}s^6+r^6s^{12}-3r^6s^6$: an exposed highest square face cannot
cancel. The permitted exponents are exactly

$$
E=\{(0,0),(2,2),(2,3),(2,4),(3,2),(3,3),(3,4),
(3,6),(4,2),(4,3),(4,4),(6,3)\}.
$$

Average the Gram under the signs of $r,s$. Successive zero coefficients
force the following root coefficients to vanish:

| Polynomial coefficient | Forced root exponent |
|---|---|
| $(4,8)$ | $(2,4)$ |
| $(4,6)$ | $(2,3)$ |
| $(8,4)$ | $(4,2)$ |
| $(6,4)$ | $(3,2)$ |
| $(6,8)$ | $(3,4)$ |
| $(8,6)$ | $(4,3)$ |
| $(8,8)$ | $(4,4)$ |
| $(4,4)$ | $(2,2)$ |

At each step the coefficient is a sum of squares of the indicated
root coefficients; every possible cross term has already vanished.
Only $(0,0),(3,3),(3,6),(6,3)$ remain. The coefficient of $r^6s^6$ is
then a sum of squared $(3,3)$ coefficients, contradicting its value $-3$.
Evaluating a state-linear certificate at $Z=1,W=X=Y=0$ gives the same
contradiction.

The [exact checker](checks/sixth_root_features.py) expands (1) and (2),
enumerates the Newton support, and verifies every forced-zero step.
The example proves that some additional sixth-root functions can be
essential. Whether the special native C1/C2 matrix identities admit a
smaller span, or admit the full $Z_6$ uniformly, remains open.
