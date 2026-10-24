"""Check 6: CHSH operator on real qubits, basis psi=(W,X,Y,Z) = |00>,|01>,|10>,|11>."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from sympy.physics.quantum import TensorProduct as kron

I2, Zp, Xp = eye(2), Matrix([[1, 0], [0, -1]]), Matrix([[0, 1], [1, 0]])
Yr = Matrix([[0, -1], [1, 0]])                      # real Pauli Y (= -i sigma_y)
ca, sa, cb, sb = (1 - u**2)/Du, 2*u/Du, (1 - v**2)/Dv, 2*v/Dv
A0 = B0 = Zp
A1, B1 = ca*Zp + sa*Xp, cb*Zp + sb*Xp
Bop = kron(A0, B0) + kron(A0, B1) + kron(A1, B0) - kron(A1, B1)
simp = lambda M: M.applyfunc(sp.cancel)
sig = sa*sb
lmb = symbols('x')

# (a)
print('(a) B^2 == 4I + 4 sa sb YxY :', simp(Bop**2 - 4*eye(4) - 4*sig*kron(Yr, Yr)) == zeros(4, 4))
print('    B^2 == 4I - 4 sa sb YxY :', simp(Bop**2 - 4*eye(4) + 4*sig*kron(Yr, Yr)) == zeros(4, 4))
# (b) charpoly of I - B/2 vs prod (x - 1 -+ sqrt(1 +- sig)) = ((x-1)^2-1-sig)((x-1)^2-1+sig)
target = ((lmb - 1)**2 - 1 - sig)*((lmb - 1)**2 - 1 + sig)
print('(b) charpoly(I-B/2) == prod(x-(1-+sqrt(1+-s))) :',
      sp.cancel((eye(4) - Bop/2).charpoly(lmb).as_expr() - expand(target)) == 0)
# (c) d-part of F = d Du Dv psi^T M psi
Fd = sp.diff(F, d)                                   # F is linear in d
M = sp.hessian(Fd, list(psi)).applyfunc(lambda e: sp.cancel(e/(2*Du*Dv)))
print('(c) tr M =', sp.simplify(M.trace()), '  tr 2(I-B/2) = 8')
print('    charpoly(M) == charpoly(I-B/2)    :', sp.cancel(M.charpoly(lmb).as_expr() - expand(target)) == 0)
target2 = ((lmb/2 - 1)**2 - 1 - sig)*((lmb/2 - 1)**2 - 1 + sig)*16
print('    charpoly(M) == charpoly(2(I-B/2)) :', sp.cancel(M.charpoly(lmb).as_expr() - expand(target2)) == 0)
Bswap = kron(A0, B0) + kron(A1, B0) + kron(A1, B1) - kron(A0, B1)   # B with A0 <-> A1
print('    M == I - Bswap/2 exactly          :', simp(M - (eye(4) - Bswap/2)) == zeros(4, 4))
# (d) every h_k is a product vector a (x) b, with A_x a = +-a and B_y b = +-b
ev = {'A0': A0, 'A1': A1, 'B0': B0, 'B1': B1}
for k in range(8):
    Hk = Matrix(2, 2, list(h[k]))                    # [[W,X],[Y,Z]]
    rank = Hk.rank(simplify=True)
    av = Hk[:, 0] if Hk[:, 0] != zeros(2, 1) else Hk[:, 1]
    bv = Hk[0, :].T if Hk[0, :] != zeros(1, 2) else Hk[1, :].T
    lab = [f'{nm}={s:+d}' for nm, Op in ev.items() for s in (1, -1)
           if simp(Op*(av if nm[0] == 'A' else bv) - s*(av if nm[0] == 'A' else bv)) == zeros(2, 1)]
    print(f'(d) {NAMES[k]:5s} rank(g g^T)={(c[k]*h[k]*h[k].T).rank(simplify=True)} Schmidt rank={rank}', lab)
