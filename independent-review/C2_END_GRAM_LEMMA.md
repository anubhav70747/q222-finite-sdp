# One three-dimensional completion at both C2 ends

At either open C2 end, completing the two state squares leaves the same
three-dimensional quadratic form. Its positive completion follows from a
single identity.

## The completion

Let \(P,Q,g,h>0\), and set
\[
 H_\tau(x,y,w)=Qx^2+Py^2-\tau(x+gw)(y-hw).
\]
Then \(H_\tau\) is positive definite exactly when
\[
 0<\tau<\tau_{\max},
 \qquad \tau_{\max}=\frac{4ghPQ}{Qg^2+Ph^2}.
\]

Indeed, completing the \(w\) square gives
\[
 H_\tau=\tau gh\left(w+\frac{hx-gy}{2gh}\right)^2
       +Ax^2-\frac{\tau}{2}xy+By^2,
 \qquad
 A=Q-\frac{\tau h}{4g},\quad B=P-\frac{\tau g}{4h}.
\]
The remaining determinant is
\[
 AB-\frac{\tau^2}{16}
 =PQ-\frac{\tau(Qg^2+Ph^2)}{4gh}.
\]
For \(0<\tau<\tau_{\max}\), both \(A\) and \(B\) are positive, and this
determinant is positive. Completing one further square proves the claim:
\[
 H_\tau=\tau gh\left(w+\frac{hx-gy}{2gh}\right)^2
       +A\left(x-\frac{\tau y}{4A}\right)^2
       +\frac{AB-\tau^2/16}{A}\,y^2.
\]
Conversely, positive definiteness requires \(\tau gh>0\) and a positive
Schur determinant, giving the same interval.

The choice \(\tau=\tau_{\max}/2\) works at every positive parameter value.
Its final determinant is \(PQ/2\). Thus no parameter search is needed.

## The two limits

The exact native-parent expansions in
[the symbolic checker](checks/gram-lift/ends/check_c2_c_limits.py) first
complete two squares \(c_X e_X^2+c_Y e_Y^2\), with \(c_X,c_Y>0\).
The remaining two squares have the form \(Qx^2+Py^2\).
The following substitutions put both ends into the completion above.

| End | \(x\) | \(y\) | \(w\) | \(g\) | \(h\) |
|---|---|---|---|---|---|
| \(c\to\infty\) | \(UW-tVZ\) | \(VW-a^2tUZ\) | \(VZ-aUZ\) | \(t\) | \(at\) |
| \(c\to0\) | \(UW-(b/a^2)VZ\) | \(VW-bUZ\) | \(VZ-aUZ\) | \(b/a^2\) | \(b/a\) |

At infinity, \(Q=q_2\) and \(P=p_2\) in the checker. At zero,
\(Q=P_0\) and \(P=Q_0\). These coefficients are strictly positive on
the respective open faces.

The subtracted product is itself a square after the monomial relations
are restored:
\[
 (x+gw)(y-hw)=
 \begin{cases}
 UV(W-atZ)^2,&c\to\infty,\\
 UV(W-(b/a)Z)^2,&c\to0.
 \end{cases}
\]
Adding \(\tau\) times this square preserves the original polynomial.
It uses the two features \(\sqrt{UV}W,\sqrt{UV}Z\).
The other block uses six ordinary monomials. The five forms
\(e_X,e_Y,x,y,w\) are independent and vanish on the contact vector,
so the six-by-six block has exactly that one-dimensional kernel.
The two-by-two block has rank one and its contact kernel.
The complete limiting Gram therefore has rank six.

For \(c\to\infty\), the face coordinates are
\[
 b=a^2t+\eta/c^2,\quad r=\rho/c^2,\qquad
 0<\rho<t,\quad \rho>t-1,\quad 0<\eta<1,\quad a^2t+\eta>1.
\]
The rescaled parent is \(c^{-2}F_c(cU,cV,W,X/c,Y/c,Z)\).
For \(c\to0\), \(a,b,t,r\) are fixed with
\[
 0<r,t<1,\qquad b+r>1,\qquad b>a^2t,
\]
and the rescaling is \(c^{-2}F_c(cU,cV,W,cX,cY,Z)\).

These are exact certificates for the two limit forms. Extending them to
finite \(c\) requires a positive completion on the additional feature
directions that appear away from the limit. Positivity on the limiting
six-dimensional range alone does not supply that extension.

