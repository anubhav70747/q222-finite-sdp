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
remains **348 pages**: its C1 calculations are grouped into seven mathematical
domains and its C2 calculations into seventeen families. The displayed
identities, theorem statements and dependencies are retained. The short
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
| Global candidates of orders 8,368 and 10,736 | Strict outer approximations, with exact separation proofs |
| Global candidate of order 13,744 | Forward inclusion; reverse inclusion remains unproved |
| C160 candidate of order 63,888 | Reverse inclusion remains unproved; full forward inclusion uses the stated extremal-coverage theorem |

The explicit-construction note derives the exact projection, gives the local
feature formulas and explains the global candidates. It also formulates an
exact algebraic search that terminates if a finite lift exists. That search
has not produced a full-set coefficient tuple.

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
