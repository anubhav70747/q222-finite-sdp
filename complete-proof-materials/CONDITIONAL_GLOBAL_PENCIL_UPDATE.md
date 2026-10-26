# Evaluated conditional global SDP

The recovered package now gives every rational coefficient of a global outer SDP for Q222. The independent analytical and exact checks support its full matrix specification. Equality with the quantum set remains conditional on H1 and H2: direct K364 admission of every physical oriented strict C1 source and every physical strict C2 source. Neither admission theorem is proved by this package.

## Delivered matrices

| Data | Exact size |
|---|---:|
| Fixed testing/state features | 364 |
| Polynomial coefficient keys | 3,250 |
| Source map T | 3,250 by 9, with 120 nonzero entries |
| Independent upper Gram coordinates | 66,430 |
| Physical relabelings | 128 |
| PSD blocks, including the event block | 897 |
| Total PSD order | 51,088 |
| Auxiliary scalar variables | 423,168 |
| Physical Bell coordinates | 8 |
| Affine equalities | 3,456 |
| Event inequalities | 16 |

The data include the full source map, every local Gram/adjoint entry, all physical action and target matrices, both closed-face source/dual bases, every affine equality and all variable/block offsets. Repeated moment blocks are specified by their complete local entries and explicit disjoint translations. Expanding those repeated entries into an additional 8.5-million-entry file is unnecessary to evaluate the matrices. The precise formulation is a block-diagonal PSD constraint together with 3,456 affine equalities; it is not an equality-free LMI of the same reported order.

For each physical action g, the target is b_g(q)=E_SG1^T Q_g(1,q). The K364 block imposes b_g=T^T y_g+r_g, A*(y_g)>=0 and r_g>=0. SG12 and SG13 use their retained-source injections, positive orthants and free-auxiliary cancellation equations. All 128 physical relabelings are present, including all eight PR sectors and the C1 reflection. No nonlinear native normalization enters the affine system.

The original specification's section 3 says there are 17 face equalities per group. The correct number is 18: seven SG12 target equations, three SG12 auxiliary cancellations, seven SG13 target equations and one SG13 auxiliary cancellation. With nine K364 equations this gives the correctly materialized 27 per group and 3,456 overall. The original file is preserved unchanged; this prose typo does not affect the data.

## Verified inclusion and remaining hypotheses

Every represented source cone is globally sound and its source dual is attained. Consequently every quantum behaviour has feasible auxiliary variables in this outer SDP. The imported original-Z6 strict witness is serialized and verified; it supplies the K364 Slater point without any universal Z6 admission assumption.

Under H1 and H2, Kbar=K364+sigma(K364) contains all strict contact sources. The verified coefficient-space density theorem then supplies the permitted singular sources. The full physical orbit absorbs sigma, while the retained event and closed-face cones cover the other source types. This proves the reverse inclusion conditional on H1/H2. The final exactness claim remains open solely because those two admission inputs are unproved for this candidate; the shorter complete proof and manuscript delivery also remain unfinished.

## Verification and preservation

Both original ZIPs and all 49 archived members pass CRC and byte/hash checks. The C1 ZIP contains one historical Python bytecode file: it is preserved inside that unchanged ZIP but omitted from the expanded public source copies. Use the ZIP for a complete replay of its original archive manifest. Every scientific source, proof and data member is published.

The supplied exact checker and generation/assembly scripts were replayed once in private copies. All 12 generated matrix JSON/JSON.gz files agree semantically with the originals. Three are byte-identical; the others differ only in terminal newline or compressed-file representation. These are semantic regeneration checks, not a claim that every regenerated file has identical bytes.

A separate verifier, without importing the supplied generator or assembler, checks all 66,430 local Gram/adjoint entries and their symmetric factors, all 128 copy partitions and offsets, all 3,456 affine rows with 34,304 nonzero entries, and all 16 events. Independent analytical review checks the source interpretation, face duals, strict feasibility, physical symmetry and conditional inclusion argument. These checks do not establish H1 or H2.

## C1 charge result

The accompanying C1 theorem concerns a fixed a=t=1, k=sqrt(x), g=sqrt(z) subfamily, with 0<x<=1, 0<z<1 and x*z^2-(1+x)^2*z+x*(1+x)>0. For each fixed such pair, sufficiently small f>0 gives strict native sources whose scalar charge demand and capacity inequalities admit gamma=0 and delta=S_Y. These charges are nonnegative, not both positive. No uniform cutoff up to the source boundary is supplied. The seven remaining gauges and two full Schur PSD conditions are unresolved, so this is not C1 admission even on that subfamily. The qualitative V01-V06 conjunction is unchanged.

## Evidence

- [Original proof and matrix specification](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819/task7/package/PROOF_AND_MATRIX_SPECIFICATION.md)
- [Original full matrix/checker package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819/task7/original-package.zip)
- [Complete indexed rational pencil](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819/task7/package/K364_FULL128_FACTORED_SPARSE_PENCIL.json.gz)
- [Independent matrix and inclusion review](materials/historical/analysis/global_pencil_analytic_verification_2026_10_26.txt)
- [Exact replay, regeneration and independent indexing](materials/historical/analysis/global_pencil_exact_verification_2026_10_26.txt)
- [Independent all-row and offset checks](materials/historical/analysis/global-pencil-exact-verification-2026-10-26/independent_indexing_results.json)
- [Original fixed-subfamily charge theorem](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819/task5/package/Q222_C1_ORIENTED_K364_ADMISSION_2026_10_26/PROOF.md)
- [Original complete C1 charge package](materials/historical/analysis/exact-lift-2026-10-23/seven-extra-high-2026-10-26/fetched-1819/task5/original-package.zip)
- [Independent small-f charge review](materials/historical/analysis/c1_small_f_charge_analytic_verification_2026_10_26.txt)
- [Targeted source-charge replay](materials/historical/analysis/global_pencil_exact_verification_2026_10_26.txt)
