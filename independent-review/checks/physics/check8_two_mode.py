"""Two-mode normal form of a two-setting Bell operator on real qubits.

Bell basis: Phi+ = (|00>+|11>)/r2, Psi- = (|01>-|10>)/r2  (parity Gamma = Y(x)Y = +1)
            Phi- = (|00>-|11>)/r2, Psi+ = (|01>+|10>)/r2  (parity -1)
With z+ = <Phi+|psi> + i<Psi-|psi> and z- = <Phi-|psi> + i<Psi+|psi>, every Bell operator
f0 + sum f_xy A_x B_y + sum f_x A_x + sum g_y B_y, with A_x = cos(al_x) Z + sin(al_x) X and
B_y = cos(be_y) Z + sin(be_y) X, has the quadratic form
  f0 (|z+|^2 + |z-|^2) + Re(W+ conj(z+)^2) + Re(W- conj(z-)^2) + 2 Re(Om z+ conj(z-)) + 2 Re(N z+ z-),
  W+ = sum f_xy e^{i(be_y - al_x)},  W- = sum f_xy e^{i(al_x + be_y)},
  Om = sum f_x e^{i al_x},  N = sum g_y e^{-i be_y}."""
import itertools
import sympy as sp
from sympy.physics.quantum import TensorProduct as kr

r2 = sp.sqrt(2)
al = sp.symbols('alpha0 alpha1', real=True); be = sp.symbols('beta0 beta1', real=True)
f0 = sp.Symbol('f0', real=True)
fxy = sp.Matrix(2, 2, lambda x, y: sp.Symbol(f'f{x}{y}', real=True))
fx = sp.symbols('fa0 fa1', real=True); gy = sp.symbols('gb0 gb1', real=True)
I2, Zp, Xp = sp.eye(2), sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, 1], [1, 0]])
obs = lambda t: sp.cos(t)*Zp + sp.sin(t)*Xp
A = [obs(a) for a in al]; B = [obs(b) for b in be]
Fop = f0*sp.eye(4) + sum((fxy[x, y]*kr(A[x], B[y]) for x in range(2) for y in range(2)), sp.zeros(4))
Fop += sum((fx[x]*kr(A[x], I2) for x in range(2)), sp.zeros(4)) + sum((gy[y]*kr(I2, B[y]) for y in range(2)), sp.zeros(4))
k00, k01, k10, k11 = [sp.Matrix([1 if i == j else 0 for i in range(4)]) for j in range(4)]
basis = [(k00 + k11)/r2, (k01 - k10)/r2, (k00 - k11)/r2, (k01 + k10)/r2]   # Phi+, Psi-, Phi-, Psi+
x1, y1, x2, y2 = sp.symbols('x1 y1 x2 y2', real=True)
psi = x1*basis[0] + y1*basis[1] + x2*basis[2] + y2*basis[3]
Q = sp.expand((psi.T*Fop*psi)[0])

zp, zm = x1 + sp.I*y1, x2 + sp.I*y2
E = lambda t: sp.exp(sp.I*t)
best = None
for s1, s2, s3, s4 in itertools.product([1, -1], repeat=4):
    Wp = sum(fxy[x, y]*E(s1*(al[x] - be[y])) for x in range(2) for y in range(2))
    Wm = sum(fxy[x, y]*E(s2*(al[x] + be[y])) for x in range(2) for y in range(2))
    Om = sum(fx[x]*E(s3*al[x]) for x in range(2)); N = sum(gy[y]*E(s4*be[y]) for y in range(2))
    cand = (f0*(zp*sp.conjugate(zp) + zm*sp.conjugate(zm)) + sp.re(Wp*sp.conjugate(zp)**2)
            + sp.re(Wm*sp.conjugate(zm)**2) + 2*sp.re(Om*zp*sp.conjugate(zm)) + 2*sp.re(N*zp*zm))
    diff = sp.simplify(sp.expand(sp.expand_complex(cand) - Q, trig=True))
    if diff == 0:
        best = (s1, s2, s3, s4); break
print('two-mode normal form holds with phase signs (W+, W-, Om, N) =', best)
# Tsirelson / TLM for correlation-only operators: each sector's form f0|z|^2 + Re(W conj(z)^2) has eigenvalues f0 +- |W|
w = sp.Symbol('w', positive=True); th = sp.Symbol('th', real=True)
Msec = sp.Matrix([[f0 + w*sp.cos(th), w*sp.sin(th)], [w*sp.sin(th), f0 - w*sp.cos(th)]])
print('sector eigenvalues:', [sp.simplify(e) for e in Msec.eigenvals()])
