# An exact W2 certificate at an unequal-angle C2 contact

The source below is an extreme ray of the full Bell polar and has an exact
sum-of-squares certificate in the thirteen-word space $W_2$. Its contact
state is maximally entangled, and the two measurement angles differ. The
certificate is given by the rational matrices below.

## The Bell slack and its contact

Let $A_0,A_1,B_0,B_1$ be self-adjoint reflections, with
$[A_x,B_y]=0$. Define

$$
H=I-\frac{A_0B_0-A_0B_1+A_1B_0+A_1B_1}{2},\qquad
4\Pi_{ab|xy}=(I+aA_x)(I+bB_y).
$$

The event order is

$$
\begin{array}{c|cccccccc}
 &\lambda&n&m&C&A&D&\alpha&\beta\\ \hline
(a,b|x,y)&(++|01)&(+-|00)&(+-|11)&(-+|10)&(-+|00)&(-+|11)&(--|01)&(+-|10).
\end{array}
$$

For a source $c=(d,\lambda,n,m,C,A,D,\alpha,\beta)$, put
$f_c=dH+4\sum_e c_e\Pi_e$. Here

$$
\lambda=\beta=0,\qquad
(d,n,m,C,A,D,\alpha)
=\left(1,\frac{43}{299},\frac{77}{1495},\frac{891}{598},
\frac{385}{299},\frac{119}{598},\frac9{598}\right).
\tag{1}
$$

In the symmetric native parametrization of the
[unequal-angle contact equations](C2_CONTACT_GEOMETRY.md), this is
$(a,r,s)=(2,1/2,1/16)$. Equivalently, the original native coordinates are
$(a,b_{\rm source},c,t_{\rm source},r)=(2,2,1/4,1/4,1/2)$.
All retained event coefficients are positive. With $z=c^2=1/16$, the
native quantities

$$
h=b_{\rm source}-a^2t_{\rm source},\qquad
g=(1+z)b_{\rm source}-(1-r+a^2zt_{\rm source}),\qquad
j=1+zr-t_{\rm source}
$$

give gaps $(h,g,j,1-zh,b_{\rm source}-a^2zr)$ equal to
$(1,25/16,25/32,15/16,15/8)$.

Use the real qubit observables

$$
A_0=B_0=Z,\qquad
A_1(u)=\frac{(1-u^2)Z+2uX}{1+u^2},\qquad
B_1(v)=\frac{(1-v^2)Z+2vX}{1+v^2}.
$$

At

$$
u_c=\frac14,\qquad v_c=\frac12,\qquad
\psi_c=\left(1,-\frac18,\frac18,1\right)^T,
\tag{2}
$$

$f_c\psi_c=0$, and both scalar angular derivatives
$\psi_c^T(\partial_u f_c)\psi_c$ and
$\psi_c^T(\partial_v f_c)\psi_c$ vanish. The coefficient matrix of this
unnormalized state satisfies $\Psi_c\Psi_c^T=(65/64)I$, so its normalized
state is maximally entangled.

## The exact universal certificate

The word space is

$$
W_2=(I,A_0,A_1,B_0,B_1,A_0A_1,A_1A_0,B_0B_1,B_1B_0,
       A_0B_0,A_0B_1,A_1B_0,A_1B_1).
$$

In this space define the ordered root column

$$
h=\begin{pmatrix}
A_1-I\\ B_1-I\\ A_0B_0-I\\ A_1B_1-I\\ B_0-A_0\\
A_0B_1-A_0\\ A_1B_0-A_0\\ A_0A_1-A_0\\ A_1A_0-A_0\\
B_0B_1-B_0\\ B_1B_0-B_0
\end{pmatrix}.
\tag{3}
$$

The last four entries need not be self-adjoint. The adjoint in the
following identity is essential. Set

$$
N_*=\begin{pmatrix}
\frac{1}{6} & \frac{17}{10} & \frac{6}{5} & \frac{6}{5} & \frac{4}{3} & \frac{17}{10} & \frac{221}{150} \\
-\frac{16}{51} & -\frac{6}{5} & -\frac{22}{85} & -\frac{22}{85} & -\frac{26}{51} & -\frac{6}{5} & -\frac{58}{75} \\
-\frac{175}{68} & -\frac{15}{4} & \frac{33}{34} & -\frac{65}{34} & \frac{65}{34} & \frac{13}{4} & -\frac{13}{4} \\
1 & 0 & 0 & 0 & 0 & 0 & 0 \\
\frac{175}{204} & \frac{5}{4} & -\frac{45}{34} & -\frac{45}{34} & -\frac{65}{102} & \frac{9}{4} & \frac{13}{12} \\
0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 \\
\end{pmatrix}
\tag{4}
$$

and

$$
S_*=\begin{pmatrix}
\frac{49146259}{1004640000} & -\frac{1935923}{80371200} & \frac{288229}{4018560} & 0 & -\frac{136427}{10046400} & 0 & -\frac{552231}{26790400} \\
-\frac{1935923}{80371200} & \frac{12689303}{345644000} & -\frac{11939}{229632} & -\frac{422431}{975936000} & 0 & 0 & -\frac{98483}{11481600} \\
\frac{288229}{4018560} & -\frac{11939}{229632} & \frac{3015449}{12629760} & 0 & \frac{6503}{114816} & -\frac{614509}{12629760} & 0 \\
0 & -\frac{422431}{975936000} & 0 & \frac{1}{16} & 0 & -\frac{1309}{229632} & -\frac{561}{12800} \\
-\frac{136427}{10046400} & 0 & \frac{6503}{114816} & 0 & \frac{5}{48} & -\frac{1}{32} & \frac{3}{52} \\
0 & 0 & -\frac{614509}{12629760} & -\frac{1309}{229632} & -\frac{1}{32} & \frac{1}{40} & 0 \\
-\frac{552231}{26790400} & -\frac{98483}{11481600} & 0 & -\frac{561}{12800} & \frac{3}{52} & 0 & \frac{11}{128} \\
\end{pmatrix}.
\tag{5}
$$

Then the following is an identity in the universal reflection algebra:

$$
\boxed{\quad f_c=(N_*^T h)^\dagger S_*(N_*^T h).\quad}
\tag{6}
$$

It uses only $A_x^2=B_y^2=I$ and opposite-party commutation. To check it,
move all Alice letters before all Bob letters, cancel adjacent repeated
letters within each party, and compare the resulting coefficients. There
are 41 reduced words in the expansion; the nine Bell coefficients equal
those of (1), and all 32 remaining coefficients are zero. This calculation
does not use a qubit representation.

For completeness, exact $LDL^T$ elimination without pivoting proves
$S_*>0$. With $L_{ii}=1$, the recursion is

$$
\delta_i=S_{ii}-\sum_{j<i}L_{ij}^2\delta_j,\qquad
L_{ki}=\frac{S_{ki}-\sum_{j<i}L_{kj}L_{ij}\delta_j}{\delta_i}
\quad(k>i).
$$

Its seven rational pivots satisfy, componentwise,

$$
(\delta_1,\ldots,\delta_7)>
\left(\frac1{25},\frac1{50},\frac1{10},\frac1{20},
      \frac1{20},\frac1{250},\frac1{1000}\right).
\tag{7}
$$

The [standalone exact checker](checks/c2_w2_unequal.py) embeds both
matrices, carries out this recurrence, verifies $L\operatorname{diag}
(\delta_i)L^T=S_*$, and checks every coefficient of (6). Its
[exact output](checks/c2_w2_unequal.json) records the rational pivots.
Thus (6) is a sum of seven operator Hermitian squares. In particular,
$f_c\ge0$ for arbitrary Hilbert-space dimension and arbitrary independent
reflection observables.

The contact reduction can also be checked directly. Let the $i$th column
of $T$ be $h_i(u_c,v_c)\psi_c$. Columns $1,2,3,5$ have determinant
$-15/136$, while $TN_*=0$. Rows $4,6,7,8,9,10,11$ of $N_*$ form $I_7$.
Therefore $N_*$ is a basis of the complete contact kernel of $T$.

## Extremality in the full Bell polar

Both deterministic assignments

$$
v_\lambda=(1,1,1,1),\qquad
v_\beta=(-1,1,-1,1)
$$

are zeros of $f_c$, in the order $(A_0,A_1,B_0,B_1)$. Consider any Bell
operator $F=\sum_{j=1}^9 b_j O_j$, with

$$
(O_1,\ldots,O_9)
=(I,A_0,A_1,B_0,B_1,A_0B_0,A_0B_1,A_1B_0,A_1B_1).
$$

The four components of $F(u_c,v_c)\psi_c$, its two scalar angular
stationarity conditions, and its evaluations at $v_\lambda,v_\beta$
form eight homogeneous linear equations in $b$. The resulting $8\times9$
matrix has a first-eight-column determinant

$$-\frac{2970}{4913}\ne0.\tag{8}$$

Here the operators are the normalized Pauli observables displayed above,
and the state in (2) is left unnormalized. The checker reconstructs this
matrix directly and verifies that the Bell coefficient vector of (1)
lies in its kernel.

Suppose $f_c=F_1+F_2$ with both $F_i$ quantum-valid Bell slacks. Each
summand must vanish at the two deterministic zeros and at the quantum
contact. Positivity at the contact gives $F_i\psi_c=0$; nonnegativity for
all nearby angles gives its two scalar stationarity conditions. Hence
both coefficient vectors obey the same rank-eight system and are
proportional to $f_c$. These proportionality constants are nonnegative,
since $f_c$ has a positive value at the uniform behavior. This proves
that (1) generates an extreme ray of the full nine-coordinate Bell
polar, rather than merely of the canonical nonnegative source cone.

Extending this certificate to a parameter neighborhood requires the
universal coefficient identities to remain solvable as the contact varies.
That additional condition remains open.
