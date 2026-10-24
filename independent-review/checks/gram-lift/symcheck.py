import sympy as sp
u, v, W, X, Y, Z = sp.symbols('u v W X Y Z', real=True)
lam, n, m, C, A, D, al, be, d = sp.symbols('lam n m C A D alpha beta d', real=True)
def F(u, v, W, X, Y, Z):
    Du, Dv = 1 + u**2, 1 + v**2
    return (4*lam*Du*(W + v*X)**2 + 4*n*Du*Dv*X**2 + 4*m*(v*W - X + u*v*Y - u*Z)**2
            + 4*C*Dv*(u*W - Y)**2 + 4*A*Du*Dv*Y**2 + (4*D + 2*d)*(u*W + u*v*X - Y - v*Z)**2
            + (4*al + 2*d)*Du*(v*Y - Z)**2 + (4*be + 2*d)*Dv*(X + u*Z)**2 - 2*d*Du*Dv*Z**2)
F0 = F(u, v, W, X, Y, Z)
print('u-flip (u,Y,Z)->(-u,-Y,-Z): diff =', sp.expand(F(-u, v, W, X, -Y, -Z) - F0))
print('v-flip (v,X,Z)->(-v,-X,-Z): diff =', sp.expand(F(u, -v, W, -X, Y, -Z) - F0))
