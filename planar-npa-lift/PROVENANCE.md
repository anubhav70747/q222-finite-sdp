# Source provenance

The article and note develop the common-tilt planar lift from the six-page
analytical extension prepared on 27 October 2026. All support signs,
boundary cases, explicit dual matrices, the nonpolyhedrality proof and
signed asymmetric branches are retained in Part I of the complete proofs.
Section F5 supplies the physical proof for all 32 projections.

The original finite-lift proof, matrices and delivery archives are
[preserved separately](https://github.com/anubhav70747/q222-finite-sdp/tree/main/extracted-lifts/common-tilt-planar).

Part II reorganizes the complete argument of Sections S1–S10 of the prior
PRL supplement into five sections. All substantive results and their
proofs are retained; repeated explanations and duplicate derivations are
merged. The original source files remain unchanged. The preceding
[23-page companion](https://github.com/anubhav70747/q222-finite-sdp/blob/3de9b29f3d09be1397124c870c31380a2a6a182d/planar-npa-lift/proofs/main.pdf)
and its
[source](https://github.com/anubhav70747/q222-finite-sdp/blob/3de9b29f3d09be1397124c870c31380a2a6a182d/planar-npa-lift/proofs/main.tex)
remain available in Git history. The current companion is 14 pages at the
same 11-point font and 23-mm margins.

| Prior argument | Current location |
|---|---|
| Full correlations, quantum models, compactness, finite NPA and dual attainment; trace bound and exact finite-level examples | S1 |
| Physical Bell matrix, excess quadratic form, exact quartic, endpoint root, cubic and higher-order expansions, detector loss and coefficient 54 | S2 |
| Independent-angle trial vector, uniform gap expansion, axis values, shifted Motzkin contradiction and equal-angle restriction | S3 |
| Uniform Taylor approximation, fixed-grid coefficient bound, polynomial sum-of-squares limit, fixed-level divergence and its quantifiers | S4 |
| Valid probabilities, arbitrary finite projector-word lists, strict inclusion and deterministic-limit mixtures | S5 |

The Letter has two pages of argument and one references page. It includes the
characteristic-coefficient test, AM–GM reduction, explicit local sum of
squares, univariate Gram construction and both physical membership inclusions.
The complete proof companion gives the actual matrices and full NPA proof.

Scientific inputs retain their citations: Hanselka–Sinn, Corollary 2.8
(including singular PSD matrix polynomials); Scheiderer, Theorem 0.1
(closed convex semialgebraic subsets of the plane); and Gigena et al.
(2025), for the symmetric tilted quantum bound. The studies of parallel
repetition and BoundingContextuality informed exposition only and are not
redistributed here.

Source hashes (SHA-256):

- finite-extension/main.tex: d5778e636fb5825e1697f16b8293a518ba8a4a49469c7900a5e36f8500ecb0e8
- prior-PRL/supplement.tex: f54cc00aa3f2b9fcad44e2377b06c491d03a1b5f0a78bde9b8aa4e7b29aac0e1

The author explicitly accepted the conceptual abstract before this
revision. ABSTRACT.txt and the article reproduce it verbatim. Its UTF-8
SHA-256, including the final newline, is
b6f96dca28746fe5f6c7c8f0ca120832819e611ae24b151fa492cca23d6c6854.

A bounded Sol high writing pass condensed Part II with a result-by-result
preservation map. The lead checked its mathematical dependencies against
the original proof. A separate Luna high reader received the Letter in
two sequential frozen portions and could reconstruct each without later
text supplying a missing inference. These are writing and proof-scope
reviews; they do not imply coauthor acceptance of the revised article.

## 28 October physical-witness revision

The previous 13-page companion is preserved at commit
d01e5e089c5b4ad73667e922eec5b62b5d82485d. Its entire Part II is unchanged
in this revision. Part I adds an explicit homogeneous state-to-membership
map, a direct proof of both inclusions, and an elementary connected-set
chord argument. Together they prove that two real pure equal-overlap
two-qubit strategies suffice for every point of the planar region, with
simultaneous witness ranks at most two. The earlier strict-anchor duality
proof is retained as an alternative.

The parallel-repetition and BoundingContextuality coauthored sources were
read again for exposition: explain the object needed, construct it, then
verify its defining properties. No isolated passage is attributed to one
coauthor and no coauthor approval is inferred. The accepted abstract is
unchanged.

Two bounded Sol high reviewers were used: one audited the science and the
new physical/chord argument; one read the whole Letter sequentially and
identified the first missing explanations. A focused return reading of
the replacement membership passage checked whether both inclusions could
now be reconstructed. This was a whole-manuscript reader followed by an
informed return pass, not a new frozen-prefix experiment.
