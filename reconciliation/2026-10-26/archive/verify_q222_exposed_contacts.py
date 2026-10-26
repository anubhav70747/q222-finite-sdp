"""Exact-rational check of the Q222 C1/C2 contact-minor mechanism.

Uses only the Python standard library.  The analytic universal nonvanishing
of the minors is EN8/EN16 in the q222-finite-sdp repository; this script
independently verifies the stationary identities and rank at two rational
seeds, not the universal claim.
"""

from fractions import Fraction as Q


class Dual:
    def __init__(self, value, gradient=None):
        self.value = Q(value)
        self.gradient = tuple(Q(0) for _ in range(5)) if gradient is None else tuple(gradient)

    @staticmethod
    def cast(x):
        return x if isinstance(x, Dual) else Dual(x)

    @staticmethod
    def coordinate(x, i):
        gradient = [Q(0)] * 5
        gradient[i] = Q(1)
        return Dual(x, gradient)

    def __add__(self, rhs):
        rhs = self.cast(rhs)
        return Dual(self.value + rhs.value, [a + b for a, b in zip(self.gradient, rhs.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-a for a in self.gradient])

    def __sub__(self, rhs):
        return self + (-self.cast(rhs))

    def __rsub__(self, lhs):
        return self.cast(lhs) - self

    def __mul__(self, rhs):
        rhs = self.cast(rhs)
        return Dual(
            self.value * rhs.value,
            [self.value * b + rhs.value * a for a, b in zip(self.gradient, rhs.gradient)],
        )

    __rmul__ = __mul__

    def __pow__(self, power):
        assert power == 2
        return self * self

    def jet(self):
        return (self.value, *self.gradient)


SOURCE_ORDER = ("d", "lambda", "n", "m", "C", "A", "D", "alpha", "beta")


def source_jets(u, v, W, X, Y):
    # Derivative order is (W,X,Y,u,v); Z is fixed to one.
    W, X, Y, u, v = (Dual.coordinate(z, k) for k, z in enumerate((W, X, Y, u, v)))
    du, dv = 1 + u**2, 1 + v**2
    event_D = u * W + u * v * X - Y - v
    event_alpha = v * Y - 1
    event_beta = X + u
    return {
        "d": 2 * event_D**2 + 2 * du * event_alpha**2
             + 2 * dv * event_beta**2 - 2 * du * dv,
        "lambda": 4 * du * (W + v * X)**2,
        "n": 4 * du * dv * X**2,
        "m": 4 * (v * W - X + u * v * Y - u)**2,
        "C": 4 * dv * (u * W - Y)**2,
        "A": 4 * du * dv * Y**2,
        "D": 4 * event_D**2,
        "alpha": 4 * du * event_alpha**2,
        "beta": 4 * dv * event_beta**2,
    }


def determinant(matrix):
    matrix = [[Q(x) for x in row] for row in matrix]
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    prod = Q(1)
    for j in range(n):
        i = next((k for k in range(j, n) if matrix[k][j]), None)
        if i is None:
            return Q(0)
        if i != j:
            matrix[j], matrix[i] = matrix[i], matrix[j]
            prod = -prod
        pivot = matrix[j][j]
        prod *= pivot
        for i in range(j + 1, n):
            z = matrix[i][j] / pivot
            for k in range(j + 1, n):
                matrix[i][k] -= z * matrix[j][k]
            matrix[i][j] = 0
    return prod


def check(name, source, contact, selected):
    jets = source_jets(*contact)
    assert len(jets) == 9
    total = [sum(source[label] * jets[label].jet()[i] for label in SOURCE_ORDER)
             for i in range(6)]
    assert total == [Q(0)] * 6, (name, total)
    matrix = [[jets[label].jet()[i] for label in selected] for i in range(6)]
    minor = determinant(matrix)
    assert minor != 0
    print(f"{name}: exact six stationary equations PASS; selected 6x6 determinant = {minor}")
    u, v, W, X, Y = contact
    w = (1 + u*u) * (1 + v*v) * (W*W + X*X + Y*Y + 1)
    h_at_contact = jets["d"].value / w
    s_at_contact = 2 * (1 - h_at_contact)
    s_at_mixture = 2 + (w / (w + Q(1, 2))) * (s_at_contact - 2)
    expected_s = Q(107, 39) if name == "C1" else Q(8, 3)
    assert s_at_mixture == expected_s
    print(f"{name}: w={w}; contact CHSH={s_at_contact}; exposing-mixture CHSH={s_at_mixture} PASS")


if __name__ == "__main__":
    C1 = dict(zip(SOURCE_ORDER, map(Q, (
        "27/2", "0", "0", "4", "2", "9", "5/4", "13/4", "9/4"
    ))))
    C2 = dict(zip(SOURCE_ORDER, map(Q, (
        "1", "0", "1/8", "3/8", "1/6", "5/8", "3/8", "1/6", "0"
    ))))
    check("C1", C1, tuple(map(Q, ("1", "1", "1", "-1/2", "1/4"))),
          ("d", "m", "A", "D", "alpha", "beta"))
    check("C2", C2, tuple(map(Q, ("1", "1", "1", "-1/4", "1/4"))),
          ("n", "m", "C", "A", "D", "alpha"))