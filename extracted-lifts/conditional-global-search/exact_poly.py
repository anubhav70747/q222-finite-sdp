"""Small exact sparse polynomial and matrix helpers (Python standard library).

Coefficients are fractions.Fraction. No floating-point arithmetic is used.
This helper is intentionally limited to polynomial expansion and substitution;
it is not a quantifier-elimination engine or a general proof assistant.
"""
from fractions import Fraction as F


class P:
    def __init__(self, value=0):
        if isinstance(value, P):
            self.d = dict(value.d)
        elif isinstance(value, dict):
            self.d = {m: F(c) for m, c in value.items() if c}
        else:
            self.d = {(): F(value)} if value else {}

    @staticmethod
    def var(name):
        return P({((name, 1),): F(1)})

    def __add__(self, other):
        other = P(other)
        d = dict(self.d)
        for m, c in other.d.items():
            d[m] = d.get(m, F(0)) + c
        return P(d)

    __radd__ = __add__

    def __neg__(self):
        return P({m: -c for m, c in self.d.items()})

    def __sub__(self, other):
        return self + -P(other)

    def __rsub__(self, other):
        return P(other) + -self

    def __mul__(self, other):
        other = P(other)
        d = {}
        for a, x in self.d.items():
            for b, y in other.d.items():
                e = dict(a)
                for n, k in b:
                    e[n] = e.get(n, 0) + k
                m = tuple(sorted(e.items()))
                d[m] = d.get(m, F(0)) + x * y
        return P(d)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * (F(1) / F(scalar))

    def __pow__(self, exponent):
        assert isinstance(exponent, int) and exponent >= 0
        ans, base = P(1), self
        while exponent:
            if exponent & 1:
                ans = ans * base
            base = base * base
            exponent >>= 1
        return ans

    def __eq__(self, other):
        return self.d == P(other).d

    def __bool__(self):
        return bool(self.d)

    def sub(self, replacements):
        ans = P()
        for m, c in self.d.items():
            term = P(c)
            for n, k in m:
                term *= P(replacements.get(n, P.var(n))) ** k
            ans += term
        return ans

    def reduce_square(self, name, replacement):
        """Use name**2 = replacement when replacement contains no name."""
        assert all(name not in dict(m) for m in P(replacement).d)
        ans = P()
        for m, c in self.d.items():
            e = dict(m)
            k = e.pop(name, 0)
            if k % 2:
                e[name] = 1
            ans += P({tuple(sorted(e.items())): c}) * P(replacement) ** (k // 2)
        return ans

    @property
    def degree(self):
        return max((sum(k for _, k in m) for m in self.d), default=0)

    def wolfram(self):
        if not self.d:
            return "0"
        terms = []
        for m, c in sorted(self.d.items()):
            factors = [str(c)] + [n if k == 1 else f"{n}^{k}" for n, k in m]
            terms.append("*".join(factors))
        return "(" + "+".join(terms) + ")"

    def json_terms(self):
        return [{"coefficient": str(c), "monomial": dict(m)}
                for m, c in sorted(self.d.items())]

    def __repr__(self):
        return self.wolfram()


def mat_add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def mat_scale(c, a):
    return [[c * x for x in row] for row in a]


def mat_mul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), P())
             for j in range(len(b[0]))] for i in range(len(a))]


def kron(a, b):
    return [[a[i][j] * b[k][l] for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def expectation(psi, operator):
    return sum((psi[i] * operator[i][j] * psi[j]
                for i in range(len(psi)) for j in range(len(psi))), P())


def rank_rational(matrix):
    a = [[F(v) for v in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        p = a[row][col]
        a[row] = [v / p for v in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                p = a[i][col]
                a[i] = [v - p * w for v, w in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row
