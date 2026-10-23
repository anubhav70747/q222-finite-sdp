# Why the earlier arc construction needs a different completion

An exact stationary zero shows why the earlier separated arc matrix fails
for the larger source family. The uniform contact certificate now covers this
zero: it lies inside the enlarged cap `|u| <= 9/256`, `|v| <= 1/16`.

Use the source dictionary in the accompanying note with

\[
B=x=2^{-20},\qquad T=2^{-40},\qquad j=1,\qquad p=H=\tfrac12.
\]

All six event weights are positive. In particular, the numerator defining
\(n\) reduces to \(\mathsf N=1/2-B/4-B^2T/4>0\).
Choose \(u=2^{-10}\), \(v=2^{-30}\), and the unnormalized real state

\[
(W,X,Y,Z)=
\bigl(1+2^{-61},\,-2^{-31},\,
2^{-10}(1-2^{-41}),\,2^{-20}(1+2^{-21})\bigr).
\]

The Bell slack vanishes exactly. This is also immediate in the rotated
coordinates: \(w=1\), \(\zeta=\sqrt{jT}\), \(t=p\sqrt y\),
\(r=-H\sqrt{yjT}\) give \(U=V=Z_*=1\) and both normals zero.
Every term of the complete parent then vanishes.

This testing point lies inside
\(2^{-24}<u^2<2^{24}\), \(|v|<2^{-25}\).
It also lies in the Alice arc used by the
[earlier regional construction](../fixed-c2-gates/Fixed_C2_Completions.pdf):
\(b=(1+x)^{-1/2}\leq1-x/4<1-2^{-26}\).
The Bob complement does not apply there.

That arc construction leaves the two terms \(nX^2\) and \(C(cE)^2\)
as separate nonnegative squares, with
\(c=v/\sqrt{1+v^2}\) and \(E=(uW-Y)/\sqrt{1+u^2}\).
At this zero both are strictly positive:

\[
X^2=(j-p)^2y>0,\qquad
(cE)^2=H^2T^2xy\,\frac{1+x}{1+y}>0.
\]

A positive Gram for the remaining terms would therefore make the slack
strictly positive, a contradiction. No choice of that construction's free
auxiliaries repairs this source. An extension must combine these event
directions with the other state directions so that every square vanishes at
the contact. The combined-normal certificate does so on the enlarged cap.

`check_annulus_obstruction.py` verifies the zero, strict inequalities and
source coefficients by exact rational arithmetic. The original argument for
infinitesimal sources places the contact below every positive real annulus
cutoff; that argument alone gives no annulus certificate for the larger
fixed numerical source box used here.
