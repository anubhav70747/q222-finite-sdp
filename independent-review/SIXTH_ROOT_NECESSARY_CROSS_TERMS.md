# Sixth-root certificates require complementary fractional cross terms

The exact strict native C2 source and its 90 integer moments are given in
[the half-power obstruction](C2_HALF_POWER_OBSTRUCTION.md). This note gives
a necessary condition on every sixth-root Gram for that source. It does
not establish or refute sixth-root admission.

Use the fixed contact normalization

$$
\widehat F(u,v,\phi)
=F_c(u_cu,v_cv,\operatorname{diag}(\psi_c)\phi),
\qquad \widehat F(1,1,\mathbf1)=0,
$$

where

$$
(u_c,v_c;\psi_c)=
\left(\frac{47}{15},\frac{611}{9000};
\frac{303}{65},-\frac{611}{2310},\frac{6621971}{455000},1\right).
$$

The testing angles satisfy $u,v\geq0$, and the four components of $\phi$
are independent real variables. Start with
$z_{ij\ell}=u^{i/6}v^{j/6}\phi_\ell$, $0\leq i,j\leq6$.
The Newton polygons of the $WW$ and $ZZ$ coefficients force every Gram
row with, respectively, $i+j<6$ or $i+j>6$ to vanish. The retained
feature vector $z$ has 154 entries. Split it into

$$
\mathcal O=\{z_{ij\ell}:i,j\in\{0,3,6\}\}\cap\operatorname{supp}z,
\qquad
\mathcal N=\operatorname{supp}z\setminus\mathcal O.
$$

Thus $|\mathcal O|=30$ and $|\mathcal N|=124$. These are the original
retained monomials, before taking differences at the contact.

**Proposition.** If $\widehat F=z^TGz$ with $G\succeq0$, then, in this
fixed monomial normalization,

$$
\operatorname{tr}G_{\mathcal N\mathcal N}
\geq\frac{\delta}{124B}>9.49\times10^{-9},
\qquad
\sum_{\{e,f\}\in\mathcal C}|G_{ef}|
\geq\frac{\delta}{2B}>5.88\times10^{-7},
\tag{1}
$$

where

$$
T=24\,000\,000\,000\,000,\qquad
B=24\,002\,021\,367\,851,
$$

$$
\delta=
\frac{31917121340091771507320448459210586965990611717744093996544823253}
{1129383645613305483559298619423160947865086012060000000000}.
$$

Here $\mathcal C$ consists of unordered distinct pairs
$e=(i,j,\ell),f=(i',j',k)\in\mathcal N$ for which both $i+i'$ and
$j+j'$ are divisible by six. Every such pair contains a complementary
fractional pair $1+5=6$ or $2+4=6$ in at least one angular index.

**Proof.** Write $Y_{pq,\ell k}$ for the integer moments in the linked
half-power proof, and define a functional on the sixth-root coefficient
space by

$$
\Lambda_T(u^{p/6}v^{q/6}\phi_\ell\phi_k)=
\begin{cases}
Y_{p/6,q/6,\ell k}+T,&p,q\in\{0,6,12\},\\
0,&\text{otherwise}.
\end{cases}
\tag{2}
$$

The added term is $T$ times evaluation at $(1,1,\mathbf1)$ on the
ordinary-polynomial target. The exact half-power pairing therefore gives
$\Lambda_T(\widehat F)=-\delta$.

Let $M_{ef}=\Lambda_T(z_ez_f)$. An old angular index is a multiple of
three, whereas at least one index of a new feature is not. Hence
$M_{\mathcal O\mathcal N}=0$. Also $M_{ee}=0$ for every
$e\in\mathcal N$. Every entry satisfies $|M_{ef}|\leq B$; the largest
entry is exactly $B$.

The old block is positive definite. To see this without removing its
contact direction, order features by state $W,X,Y,Z$ and then
lexicographically by $(i,j)$. In each parity block use its first feature
as pivot. The congruence to the pivot and its differences has the form

$$
\begin{pmatrix}a_b+T&b_b^T\\b_b&Q_b\end{pmatrix},
\tag{3}
$$

where $Q_b$ is exactly the positive definite compressed moment matrix in
the half-power proof. Explicitly, write $\Lambda_0$ for (2) with $T=0$;
then $a_b=\Lambda_0(z_p^2)$,
$(b_b)_e=\Lambda_0(z_p(z_e-z_p))$ and
$(Q_b)_{ef}=\Lambda_0((z_e-z_p)(z_f-z_p))$.
Direct rational Schur complements give

| Parity | Old order | Compressed order | $a_b+T-b_b^TQ_b^{-1}b_b$ greater than |
|---|---:|---:|---:|
| $(0,0)$ | 14 | 13 | $400\,000\,000\,000$ |
| $(0,1)$ | 6 | 5 | $15\,000\,000\,000\,000$ |
| $(1,0)$ | 6 | 5 | $19\,000\,000\,000\,000$ |
| $(1,1)$ | 4 | 3 | $23\,000\,000\,000\,000$ |

Thus $M_{\mathcal O\mathcal O}\succ0$. All entries in (3) are specified
by (2) and the displayed pivot/difference prescription; the
[exact checker](checks/sixth_root_necessary_cross_terms.py) verifies these
rational inequalities from the 90 moments.

Applying $\Lambda_T$ to the Gram identity gives

$$
-\delta=
\operatorname{tr}(M_{\mathcal O\mathcal O}G_{\mathcal O\mathcal O})
+\operatorname{tr}(M_{\mathcal N\mathcal N}G_{\mathcal N\mathcal N}).
\tag{4}
$$

The first trace is nonnegative. The entry bound gives
$\|M_{\mathcal N\mathcal N}\|_{\mathrm{op}}\leq124B$, so positivity of
$G_{\mathcal N\mathcal N}$ implies the first inequality in (1).
The new block has zero diagonal and is supported on $\mathcal C$.
Consequently the absolute value of its trace pairing is at most
$2B\sum_{\mathcal C}|G_{ef}|$, proving the second inequality. The factor
two comes from the two symmetric off-diagonal entries. This argument
also applies before parity averaging, since all cross-parity entries
of $M$ vanish. $\square$

The negative pairing persists in an open neighborhood within the strict
native parameter family under the moving contact normalization. The
same fixed moment matrices therefore give positive bounds there, with
possibly decreased constants.

In particular, appending fractional squares while setting every
complementary fractional correlation to zero cannot represent this
source. The functional in (2) is **not** a separator for the sixth-root
cone: its new block is indefinite. The proposition supplies neither an
admitting sixth-root Gram nor a proof that none exists.
