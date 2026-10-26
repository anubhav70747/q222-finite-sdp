ONE SOURCE WEIGHT COMPLETES BOTH R1 CONTACT NORMALS

New proved scope
---------------
The exact original R21 RSR source dictionary with
  0 < B, x=jB, T <= 2^-10;
  0 < p < min(j,1); 0 < H < 1; n >= 0,
with no bound on j or jT, is admitted on the fixed testing cap
  |u| <= 2^-12, |v| <= 2^-25, all four real state coordinates.
The theorem includes literal n=0, every state/angle sign, testing zeros,
the positive-event cone E, and all finite homogeneous coefficient limits.

The new cancellation is phi=p^2/(1-p+p^2). It joins direct and harmonic
purchases before singular source factors are bounded. Both state normals
remain independent. A source-only 13-row core, three fixed matrix gates,
and nine fixed scalar monomials give the displayed complete certificate.

This is NOT global canonical coverage or a full Q222 lift. The first next
missing extraction is a fixed family for the actual R1 parent on a
complementary testing annulus, including n=0. The critical outer completion
uses n0>0 while n/n0=N*(j+1)/(p*P) can be zero. Other C1/C2 regions and
projective tests outside this cap also remain outside this theorem.
The original companion and the established full-polar assembly are unchanged.

Read
----
proof/all_r1_cap.pdf: 10-page proof, with 4 pages of main argument and exact
source, budget, and finite-function appendices.
proof/all_r1_cap.tex plus its three inputs: complete editable LaTeX.

Data
----
data/fixed_cap_family.json: the two actual algebraic families and all nine
  affine target columns, with explicit radical weights and monomial bases.
data/cap_affine_manifest.json: exact equation and block indexing.
data/cap_gram.coo.tsv: equation, auxiliary, rational coefficient.
data/cap_source.coo.tsv: equation, source, rational coefficient.
data/cap_pencil.coo.tsv: auxiliary, block, row, column.
All indices are zero based. Off-diagonal coefficient-equation records
already include the factor two. Pencil records are actual matrix entries.
The pencil is diag(c0,...,c8,Q0,...,Q50)>=0 plus the coefficient equations.
Total order299, 2370 free Gram entries, 554 coefficient equations,
2852 Gram records and132 source records. These counts assert no ranks.

The functions divide by sqrt(1+v^2). Their contraction equals F_h/4.
Do not divide the Gram by four again. The signed Q row is
  (v*zeta-T*u*w)/rhoQ = -sqrt(qc)*Qstar.
In complete_budget_checks.json, formal core parameter mp means mu/p;
the restored physical product m*p is mu*p/(1-p). Its meaning is recorded.

data/complete_budget_checks.json: actual 13-row coefficient records,
  source definitions, rational reserve bounds, and literal gate budgets.
data/r1_sos_roots.json: primitive rational SOS matrices with exact roots.
data/cap_closure_matrices.json: actual rational positive moment matrices
  for the fixed320-point separating grid, and the source trace functionals.

Reproduce
---------
Install the single requirement, then from any directory run:
  python /absolute/path/to/Q222_All_R1_Contact_Completion/reproduce.py
Add --pdf to compile the LaTeX twice with pdflatex.
There is no SDP solver, numerical search, or external source download.
The original R21 RFN8/9/12/21/22 parent identities are established inputs.
The new full-parent recombination is checked with both normals free.
The independent support checker keeps source scaling symbols independent.

Recovery and preservation
-------------------------
Every new proof/generator source file is shorter than18000 characters.
recovery_parts/ contains verbatim source bundles below17500 characters,
with exact file lengths and SHA-256 checksums. Data regenerate from these
sources; the bulk JSON/COO arrays need not be copied into chat.
Use recover_parts.py with the part text files to restore the sources.

The previous fixed C2 contact ZIP remains untouched. Its verified SHA-256:
c22fae5078293f1f60d8a38b55ab845740d90cad3461637663d0519fae99d522

New result only. No prior weak or strong source-affine array was regenerated.
