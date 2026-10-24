"""Parity Gamma = Y(x)Y for real qubits: commutes with A(x)B, anticommutes with marginals.
CHSH spectrum per parity sector, soft modes at compatibility, and gap to hard modes."""
import sympy as sp
from sympy.physics.quantum import TensorProduct as kr
ca, sa, cb, sb, x = sp.symbols('c_a s_a c_b s_b x', real=True)
I2, Zp, Xp, Yr = sp.eye(2), sp.Matrix([[1,0],[0,-1]]), sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-1],[1,0]])
A0, A1 = Zp, ca*Zp + sa*Xp
B0, B1 = Zp, cb*Zp + sb*Xp
G = kr(Yr, Yr)
on_circle = {ca**2: 1 - sa**2, cb**2: 1 - sb**2}
z = lambda M: M.applyfunc(lambda e: sp.simplify(sp.expand(e).subs(on_circle))) == sp.zeros(*M.shape)
print('Gamma^2 = I:', z(G*G - sp.eye(4)))
for nm, Op in [('A0B0',kr(A0,B0)),('A0B1',kr(A0,B1)),('A1B0',kr(A1,B0)),('A1B1',kr(A1,B1))]:
    print(f'[Gamma,{nm}] = 0:', z(G*Op - Op*G))
for nm, Op in [('A0',kr(A0,I2)),('A1',kr(A1,I2)),('B0',kr(I2,B0)),('B1',kr(I2,B1))]:
    print(f'{{Gamma,{nm}}} = 0:', z(G*Op + Op*G))
S = kr(A0,B0) - kr(A0,B1) + kr(A1,B0) + kr(A1,B1)          # the paper's CHSH form
sig = sa*sb
print('S^2 = 4I + 4 sin_a sin_b Gamma:', z(S*S - 4*sp.eye(4) - 4*sig*G))
# eigenvectors of Gamma, then I - S/2 restricted to each sector
Vp = sp.Matrix.hstack(*[v for val, mult, vs in G.eigenvects() if val == 1 for v in vs])
Vm = sp.Matrix.hstack(*[v for val, mult, vs in G.eigenvects() if val == -1 for v in vs])
H = sp.eye(4) - S/2
for nm, V in [('+', Vp), ('-', Vm)]:
    Hs = (V.T*H*V).applyfunc(sp.simplify)
    Ns = (V.T*V)
    cp = sp.factor(sp.expand((Hs - x*Ns).det()).subs(on_circle))
    print(f'sector Gamma={nm}1: char poly of I-S/2 (generalized) ->', cp)
print('cross block V+^T (I-S/2) V- = 0:', z(Vp.T*H*Vm))
