# Half-powers do not cover the native C2 family

An exact strict native C2 source lies outside the fixed half-power cone

$$
L_2=\{c:F_c=Z_{1/2}^TGZ_{1/2},\ G\succeq0\},\qquad
Z_{1/2}=(u^{i/2}v^{j/2}\psi_\ell)_{0\leq i,j\leq2,\,\ell=W,X,Y,Z}.
$$

The proof reduces every possible Gram to four contact blocks and gives
an integer functional positive on those blocks but negative on the
complete parent. It excludes every Gram allocation in the 36 features.

## A strict rational contact

Use the complete SG2 parent with $d=1$, $\lambda=\beta=0$ and choose

$$
(a,b,c,t,r)=
\left(\frac{13}{600},\frac{101}{1000},\frac{47}{15},
\frac{27}{35},\frac{13}{154}\right).
\tag{1}
$$

The native contact is

$$
(u_c,v_c;\psi_c)=
\left(\frac{47}{15},\frac{611}{9000};
\frac{303}{65},-\frac{611}{2310},\frac{6621971}{455000},1\right).
$$

Here are the exact stationary coefficients, specified without long
rational numerators. Put $z=c^2$, $s=1+zr$,
$h=b-a^2t$, $g=(1+z)b-(1-r+a^2zt)$ and $j=1+zr-t$, and define

$$
\begin{aligned}
\Delta={}&1+(1+z)b+z(1+z)b^2+2z(1+z)br+z^2r^2+2zr\\
&+a^2z\bigl((1+z)b+1+zr(2-r)\bigr),\\
N={}&(1+z)(b+r)(b+a^2(1-t))-a^2(1+a^2z)(1-r)t,\\
\Xi={}&(b+r)s(s+a^2zt).
\end{aligned}
$$

Solving the four state-kernel equations and two angular stationary
equations gives the native source

$$
\begin{aligned}
n&=\frac{\Xi-r\Delta}{2r\Delta},&
m&=\frac{(1+a^2z)(1-r)s(s+a^2zt)}{2g\Delta},\\
C&=\frac{sN}{2a^2t\Delta},&
A&=\frac{a^2\Xi}{2h\Delta},\\
D&=\frac{(1+z)(b+r)(b-a^2zr)(s+a^2zt)}{2a^2j\Delta},&
\alpha&=\frac{z^2(b+r)N}{2(1-zh)\Delta}.
\end{aligned}
\tag{2}
$$

All six event coefficients are strictly positive. In particular the
native gaps are

$$
\begin{gathered}
h=\frac{140893}{1400000},\quad
g=\frac{601024007}{3465000000},\quad
j=\frac{36637}{34650},\\
1-zh=\frac{3767363}{315000000},\qquad
b-a^2zr=\frac{1255020827}{12474000000}.
\end{gathered}
$$

Substitution in [SG2](../latex/chapters/04_geometry.tex) gives
$F_c=\nabla_{u,v,\psi}F_c=0$ at the displayed contact. The
stationary-contact theorem SG7 therefore proves global validity.
The [exact checker](checks/c2_half_power_obstruction.py) verifies all
these identities and signs from the full parent. Its
[output](checks/c2_half_power_obstruction.json) also gives every source
coefficient as a rational number.

## Every half-power Gram has the same forced face

Set $u=x^2$, $v=y^2$. Any polynomial SOS of $F_c(x^2,y^2,\psi)$
has roots linear in the state and of degree at most two separately
in $x,y$: highest real squares cannot cancel. Thus the 36 features
exhaust all polynomial SOS possibilities after this substitution.
Sign averaging in $x,y$ gives parity blocks of orders $(16,8,8,4)$.

The diagonal state coefficients in SG2 are

$$
\begin{aligned}
[P_c]_{WW}&=4mv^2+(4C+4D+2d)u^2+4Cu^2v^2,\\
[P_c]_{ZZ}&=4\alpha+(4\alpha+2d+4m)u^2+4Dv^2.
\end{aligned}
$$

Each root's monomial support lies in half the Newton polygon of its
square sum. This forces $i+j\geq2$ on $x^iy^jW$ and $i+j\leq2$ on
$x^iy^jZ$. Equivalently, the vanishing lowest coefficients remove
$W,xW,yW$, and the vanishing highest coefficients remove
$x^2y^2Z,x^2yZ,xy^2Z$. The remaining parity orders are $(14,6,6,4)$.

Now rescale by source constants:

$$
u=u_cU,\qquad v=v_cV,\qquad
\psi_\ell=(\psi_c)_\ell\phi_\ell.
$$

This invertible change preserves half-power membership and places the
contact at $U=V=1$, $\phi=(1,1,1,1)$. Every square vanishes there.
For each parity block choose its first retained feature $z_p$ and use
$z_e-z_p$ for every other feature. These differences span the full
contact annihilator. Every feasible Gram therefore compresses to
orders $(13,5,5,3)$, without losing a possible certificate.

## An integer functional separates the parent

For each angular pair $(p,q)$, the following row lists the moments of
the ten ordinary state monomials
$(WW,WX,WY,WZ,XX,XY,XZ,YY,YZ,ZZ)$:

```text
(0,0) 161019792 0 0 0 -121999126 -51597381 -92880749 198288996 63902153 19654676
(0,1) 187856424 -46850180 181090964 93807406 104512504 -51884026 -170947993 104746412 60911824 56523337
(0,2) 198410875 -202408890 134972214 153970197 2021367851 -131764329 -929567199 102173701 98542620 600959153
(1,0) 187856424 -20723372 228808703 94672500 -8876925 -32301218 -55239329 117094417 62928564 26333033
(1,1) 211714865 -20971443 135517649 91673664 -66738426 -13877461 -150724699 104061701 67646572 392762321
(1,2) 165685625 -100092162 132934774 129085407 872170933 -110726908 -814152246 108962482 438514151 328748742
(2,0) 259329330 -1150608 147851759 93736392 966590370 -12339807 -94532325 112427108 80542993 535327210
(2,1) 166289097 16862629 134869776 98549977 120025352 -54544777 -527190917 121915289 582567055 328748742
(2,2) 163696817 -79616153 139866971 0 703323270 -491353386 0 629878509 0 493123113
```

Denote these integers by $Y_{pq,\ell k}$. For features
$e=(i,j,\ell)$, $f=(i',j',k)$ in one parity block set

$$
\mathfrak m(e,f)=Y_{(i+i')/2,(j+j')/2,\min(\ell,k),\max(\ell,k)}.
$$

The functional on products of contact roots has matrix

$$
Q_{ef}=\mathfrak m(e,f)-\mathfrak m(e,p)
-\mathfrak m(p,f)+\mathfrak m(p,p).
\tag{3}
$$

Fix the feature order: state $W,X,Y,Z$ first, then lexicographic $(i,j)$;
use the first retained feature as pivot. In parity order
$(0,0),(0,1),(1,0),(1,1)$, exact rational $LDL^T$ factorization gives
the following lower bounds on **every** diagonal pivot:

| Block order | Lower bound |
|---|---:|
| 13 | 89 |
| 5 | 100 |
| 5 | 360 |
| 3 | 2,000,000 |

All four matrices are therefore positive definite. The checker prints
all pivots and verifies the exact matrix reconstruction. Different
feature orders change the pivots, but preserve positive definiteness.

Write $b_{pq,\ell k}$ for the ordinary polynomial coefficients of
$\widetilde F(U,V,\phi)=F_c(u_cU,v_cV,\operatorname{diag}(\psi_c)\phi)$.
Direct contraction of the entire parent with the displayed integers gives

$$
\langle Y,b\rangle=
-\frac{31917121340091771507320448459210586965990611717744093996544823253}
{1129383645613305483559298619423160947865086012060000000000}<0.
\tag{4}
$$

If a half-power SOS existed, its compressed positive semidefinite
Grams $G_j$ would satisfy
$\langle Y,b\rangle=\sum_j\operatorname{tr}(Q_jG_j)\geq0$,
contradicting (4). Off-diagonal products require no extra factor in
(3): the matrix trace supplies the two symmetric contributions.
This proves $F_c\notin L_2$.

The functional is positive on the forced contact face. Positivity on
unrestricted, uncompressed Grams is unnecessary because every Gram for
this parent must lie in that face. The same normalized matrices give
exclusion throughout an open neighborhood in the strict native C2
parameter domain: all source gaps stay positive, and (4) stays negative
by continuity. The normalization varies with the native contact.

Thus half-powers cannot provide uniform native C2 admission. Sixth-root
features include this smaller space and allow additional cubic-circuit
terms; their universal admission remains an open implication in the
[full lift construction](EXACT_LIFT_ROUTE.md).
