# Why the native C2 parent needs more than polynomial squares

The full native C2 state form has a simple necessary condition for an
ordinary polynomial sum of squares. The condition fails generically,
including at every unequal-angle maximally entangled contact. It explains
why the positive contact-scalar Gram alone cannot complete the state form.

## Four contacts force an incompatible coefficient ratio

Use the complete cleared parent $F_c(u,v,\psi)=\psi^TP_c(u,v)\psi$
from [SG2](../latex/chapters/04_geometry.tex). Set $d=1$ and
$\lambda=\beta=0$. The native parameters satisfy
$a,b,c,t>0$, $0<r<1$, $h=b-a^2t>0$, and the remaining native gaps.
Allow $n=0$. The positive contact is

$$
u_c=c,\qquad v_c=ac,\qquad
\psi_c=(W_c,X_c,Y_c,Z_c)
=\left(\frac ba,-cr,\frac{ch}{a},1\right).
$$

**Claim.** If $F_c$ is a sum of squares of real polynomials in
$(u,v,W,X,Y,Z)$, then $h=a^2r$.

The parent has degree at most two in each testing angle and is homogeneous
of degree two in the state. Each square root must therefore be of the form

$$
R_j=\sum_{\ell\in\{W,X,Y,Z\}}
(A_{j\ell}+B_{j\ell}u+C_{j\ell}v+D_{j\ell}uv)\psi_\ell.
$$

Indeed, highest-degree coefficients in a sum of squares cannot cancel.
The constant-in-state part also vanishes, since $F_c(u,v,0)=0$.
The two omitted events make the coefficients of $W^2$ and
$u^2v^2Z^2$ zero in SG2. Hence $A_{jW}=D_{jZ}=0$ for every root.

The sign symmetries

$$
(u,Y,Z)\mapsto(-u,-Y,-Z),\qquad
(v,X,Z)\mapsto(-v,-X,-Z)
$$

give four zeros of $F_c$. Every root vanishes at all four. Averaging
these four root equations retains only $W,vX,uY,uvZ$; the first and
last are absent. Thus

$$
v_cX_c C_{jX}+u_cY_c B_{jY}=0,
\qquad B_{jY}=\frac{a^2r}{h}C_{jX}.
$$

SG2 gives $[u^2Y^2]F_c=4A$ and $[v^2X^2]F_c=4n+2d$.
Comparing these coefficients in the proposed squares yields

$$
\frac{4A}{4n+2d}=\left(\frac{a^2r}{h}\right)^2.
\tag{1}
$$

The native stationary coefficients give a different relation. With
$\Xi=(b+r)(1+c^2r)(1+c^2r+a^2c^2t)>0$ and the native
denominator $\Delta>0$, they are

$$
n=\frac{\Xi-r\Delta}{2r\Delta},\qquad
A=\frac{a^2\Xi}{2h\Delta}.
$$

Since $d=1$, these imply

$$
\frac{4A}{4n+2d}=\frac{a^2r}{h}.
\tag{2}
$$

The positive ratio must equal its square, so it is one. This proves
the claim. No step divides by $n$, and the proof includes its zero face.

For example, $(a,b,c,t,r)=(1,5/4,1,3/4,1/4)$ gives

$$
(d,n,m,C,A,D,\alpha)
=\left(1,\frac7{66},\frac5{33},\frac5{22},
\frac{10}{33},\frac{16}{33},\frac9{22}\right).
$$

Here $h=1/2$ and $a^2r=1/4$: the coefficient ratio is $1/2$,
whereas the contact equations would force $1/4$ in (1).
The exceptional hypersurface $h=a^2r$ is a necessary condition;
sufficiency remains open.

For maximally entangled native contacts, $b=a$ and $t=(1-r)/a$,
so $h=ar$. The condition becomes $a=1$. Every unequal-angle
maximally entangled native parent therefore fails ordinary polynomial SOS.
This includes the source with the [exact W2 certificate](C2_W2_UNEQUAL_CERTIFICATE.md).
There is no inconsistency: if $D_u=1+u^2$ and $D_v=1+v^2$, every
$W_2$ word has common denominator $D_uD_v$. Clearing a $W_2$ SOS
certifies $D_uD_vF_c$, while the obstruction concerns $F_c$ itself.

## What remains in the coupling construction

The [contact reduction](CONTACT_SCALAR_GRAM.md) writes the state matrix as
$K=\left(\begin{smallmatrix}H&\zeta\\\zeta^T&\gamma\end{smallmatrix}\right)$,
with $\gamma=q^TGq$, $q=(s,\tau,s\tau)^T$, and
$G\succ0$ at every finite native contact. Here
$s=u-u_c$ and $\tau=v-v_c$.
For any polynomial $V$ satisfying $\zeta=Vq$, completing the square
leaves $H-VG^{-1}V^T$. If this remainder were an ordinary polynomial
matrix SOS, the full $F_c$ would be polynomial SOS as well.
The claim excludes that possibility whenever $h\ne a^2r$.

The stronger condition $VG^{-1}V^T\preceq H$ at every real test also
bounds the testing degree of $V$. Since $G$ is
constant in testing and positive definite, growth at infinity forces
$V$ to have degree at most one in each of $u,v$. For each row, the map

$$
(V_1,V_2,V_3)\longmapsto sV_1+\tau V_2+s\tau V_3
$$

has twelve bilinear coefficients, rank eight, and the four-dimensional
kernel

$$
(-\tau,s,0),\quad(-\tau,0,1),\quad
(-s\tau,0,s),\quad(0,-s\tau,\tau).
$$

These vectors are independent and each maps to zero; the image contains
all eight nonconstant monomials of bidegree at most $(2,2)$.
The three coupling rows consequently leave twelve source-side parameters.
Finding them with the required positivity, and certifying the remainder
in a fixed fractional-power space, remain open.
