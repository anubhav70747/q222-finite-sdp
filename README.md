# Finite semidefinite representation of the simplest quantum correlation set

The argument begins with a supporting Bell inequality. Its state matrix has
three positive transverse directions, and its determinant reduces to eight
angular monomials. At a stationary contact those monomials balance. Elementary
square identities certify the angular part; polynomial contact relations and
complementary source weights control the state couplings in singular limits.

- [Article](pdf/Q222_Paper.pdf): the mathematical mechanism and global argument.
- [Complete proofs](pdf/Q222_Complete_Proofs.pdf): all source domains, singular cases and complementary measurement charts.
- [Explicit representations and local certificates](pdf/Explicit_Lifts.pdf): exact formulas and the status of the recovered matrix constructions.
- [Uniform C1 certificate](pdf/Uniform_C1_Certificate.pdf): the self-contained local unification lemma.
- [Overleaf sources](Q222_Overleaf.zip).
- [Scientific matrix data and exact checkers](extracted-lifts/).

The main article is **7 pages**. The complete technical companion
is **343 pages**: its C1 calculations are grouped into seven mathematical
domains and its C2 calculations into seventeen families. One full-state completion now treats the large-$L$ C1 source for all Alice
angles, replacing three separate proofs while retaining their complete source
and testing domain. A second two-page argument now treats every testing
angle for the unit-$L$ source with unbounded $\eta$, replacing its two
separate testing proofs. The short
article uses the companion's source estimates in the canonical-certificate
proposition.

## Which matrices are explicit?

The qualitative existence argument and an evaluated matrix pencil are
different results. **An exact eight-coordinate pencil has not been extracted.**
The recovered files provide:

| Construction | What is established |
|---|---|
| Tilted-CHSH projection | Exact order-4 pencil without auxiliaries; alternative order-5 pencil with one auxiliary; both inclusions |
| Regular C2 contact box | Explicit local certificate: 785 order-2 blocks, 1,570 listed functions and 237 coefficient equations |
| C1D source box | Explicit local library: 82,928 separated squares and 3,844 testing spanning functions |
| Bounded-z C1 source box | Exact order-669 pencil for the specified local certificate cone; 120 blocks and 2,694 auxiliaries |
| Reciprocal Alice angle, every positive source ratio | Exact order-62 local source-cone pencil; 19 blocks and 97 auxiliaries |
| All Alice angles, bounded-z source box | Exact order-731 intersection of the two specified local certificate cones; 139 blocks and 2,791 auxiliaries |
| Canonical Bell validity at every angle | Exact order-49 pencil with a quartic source map; 24 order-two blocks, one scalar block and 48 auxiliaries |
| Determinant-weighted full state polynomial | Exact order-133 fixed square certificate; 104 functions and 27 blocks |
| Global candidates of orders 8,368 and 10,736 | Strict outer approximations, with exact separation proofs |
| Global candidate of order 13,744 | Forward inclusion; reverse inclusion remains unproved |
| C160 candidate of order 63,888 | Reverse inclusion remains unproved; full forward inclusion uses the stated extremal-coverage theorem |

The explicit-construction note derives the exact projection, gives the local
feature formulas and explains the global candidates. It also formulates an
exact algebraic search that terminates if a finite lift exists. That search
has not produced a full-set coefficient tuple.

The [stationary-contact construction](extracted-lifts/contact-selector/candidate_adequacy.pdf)
now supplies an exact certificate for the support that excluded the two
smaller candidates. A rational Gram selector extends it to an explicit
open five-parameter family in the order-13,744 construction. The restored
inequalities hold for all quantum realizations. This proves a local part
of the remaining converse; it does not cover every contact.

The [fixed testing-space construction](extracted-lifts/fixed-features/fixed_box_c1.pdf)
gives a complete local C1 certificate on an explicit source and testing
rectangle, including an unbounded ratio of source scales. Its two angular
cases use 220 and 157 labelled functions; their union has at most 377.
The 14-page derivation includes the complete state polynomial, zero event
weights, testing axes and all state signs. The accompanying JSON files
contain the functions, fixed directions and symbolic coefficient matrices.
This is a local construction; the remaining domains must still be covered
to extract an exact global pencil.

The [evaluated local pencil](extracted-lifts/local-order669/proof/Fixed_Family_Extraction.pdf)
uses fixed algebraic functions to represent the corresponding source cone
with rational matrices. Its 17-page derivation covers the larger source box
`0<x,sigma<=2^-30`, `0<z<=1`, under the printed gap conditions, and an
independent testing rectangle. The matrix arrays and a standard-library
checker for all 682 coefficient equations are included. This supplies an
actual local pencil while leaving the other global domains to be covered.

The [reciprocal-angle proof](extracted-lifts/inverse-alice/inverse_alice_all_z_certificate.pdf)
uses a source rescaling to turn the leading state coupling into a vector of
norm at most one. The resulting completion works for every positive source
ratio and gives 51 fixed testing functions. Its seven-page derivation includes
the exact source bounds and both directions of the order-62 pencil.
Combined with the inner-angle construction, it covers every projective Alice
angle for `0<z<=1`, `0<x,sigma<=2^-78/(1+z)` and `|w|<=2^-36`, under the
printed native gaps. The exact intersection pencil has total order 731.
The complete companion uses the same cancellation to replace its former
three C14.1 proofs with one theorem, including the remaining Bob-angle join.

The [determinant certificate](extracted-lifts/global-weighted/proof/global_weighted_certificate.pdf)
gives a seven-page argument for all nonnegative canonical source coefficients
and every measurement angle. A polynomial adjugate identity completes the
state quadratic, and twelve elementary circuits express its angular
determinant as fixed squares. The resulting order-133 pencil decides validity
of each specified canonical inequality exactly. Its input is quartic in the
inequality coefficients, so it does not supply an affine lift of the full
quantum set. The complete fixed functions, product equations, coefficient
arrays and reproduction scripts are included.

The [three-page angular-cone proof](extracted-lifts/global-weighted/angular-order49/angular_cone.pdf)
reduces the same canonical validity test to twelve balanced circuits. Each
circuit has an elementary geometric-mean bound represented by two-by-two
positive-semidefinite matrices. The resulting order-49 pencil has 48
auxiliaries and 109 rational coefficient records, all supplied explicitly.
It includes zero coefficients and minima approached only at infinity.
This pencil is affine in the nine Laurent coefficients; their dependence
on the Bell coefficients remains quartic. The order-133 construction
additionally supplies the fixed square factors for the full state polynomial.

The [bounded source-interval theorem](extracted-lifts/global-weighted/bridges/bounded_unit_L/bounded_unit_L_all_tests.pdf)
extends a fixed testing family to every pair of measurement angles when
the native source ratio belongs to a specified positive interval. Positive
event weights control the state uniformly away from the singular corner;
two endpoint completions and the corner certificate complete the cover.
The quantitative six-page proof includes exact functions, source bounds and
independent symbolic checks. Its interval-dependent grid has not been
enumerated into matrix arrays.

The [tensor-rank note](extracted-lifts/fixed-features/uniform_feature_route.pdf)
explains why a bound on each square's tensor rank differs from a common
finite testing space, and derives the final matrix construction once such
a global space is supplied.

## Build and reproduce

The editable manuscript sources are in `latex/`. With pdfLaTeX, compile each
entry twice:

```text
main.tex
Supplement.tex
Explicit_Lifts.tex
Uniform_C1_Certificate.tex
```

The bibliography is included. `main.tex` builds independently; its companion
locators refer to the accompanying version of `Supplement.tex`.

The data README gives the commands and dependencies for reproducing the exact
rational checks. The projection, regular-C2 and outer-data checks use only the
Python standard library; the C1D source and square-library checks also use
SymPy. These checks verify their printed identities and stated local domains.
They do not establish either unproved reverse inclusion.

Earlier versions remain in the Git history. The repository contains the
scientific sources, formulas, matrix arrays and reproduction code.
