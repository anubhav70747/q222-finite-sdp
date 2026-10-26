# A uniform certificate near compatible measurements

Two positive coefficients can vanish at opposite ends of the source interval.
Completing the squares separately introduces singular denominators. The weight
`phi = p^2 / (1 - p + p^2)` combines the decompositions first. It cancels both
denominators, leaving a denominator bounded below by `3/4`. Keeping the entire
first-normal coupling preserves the cancellation that makes this possible.

The resulting certificate covers the stationary source family with
`0 < B, jB, T <= 2^-10`, `0 < p < min(j,1)`, `0 < H < 1` and `n >= 0`.
There is no bound on `j` or `jT`. The testing region is
`|u| <= 2^-12`, `|v| <= 2^-25`, with every real state and both state signs.
The source parameters and testing variables are independent.

The proof retains both state directions, the literal `n = 0` case,
nonnegative event additions and finite limits of conic combinations of the
displayed sources. Fixed testing functions turn the square identities into
an affine semidefinite system. Source-dependent formulas choose its Gram
matrices; the constraints themselves are linear in the nine source
coefficients and the free Gram entries.

The complementary testing regions and the other canonical source families
remain necessary for an exact lift of the full quantum correlation set.
The earlier native-source construction and its regional matrices remain
available in [contact-completion](../contact-completion/).
