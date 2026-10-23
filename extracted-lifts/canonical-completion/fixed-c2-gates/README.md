# Fixed square decompositions away from a singular contact

[The four-page proof](Fixed_C2_Completions.pdf) derives two small affine
Gram matrices directly from the full canonical Bell slack. A cancellation
between the state normals leaves three fixed testing radicals.

- The Alice arc 0 <= b <= 1-2^-26 uses a matrix of order 20 and three
  free auxiliaries. Its 12 rational coefficient matrices have 36
  upper-triangular records.
- The complementary Bob region |c| >= 2^-26 uses a matrix of order 10,
  nine rational coefficient matrices and 15 records, with no auxiliaries.

The proof gives the source admission conditions and verifies them uniformly
on the displayed native C2 envelope, with no lower bound on its positive
openings. Both constructions include their projective testing endpoints,
all state signs and literal n=0. The nine source-sign inequalities are
separate scalar PSD blocks, giving total orders 29 and 19.

The exact native contact lies in the omitted cap b>1-2^-26, |c|<2^-26.
There the parent is zero while a retained square is strictly positive.
A different completion is required on that cap. These matrices therefore
do not constitute a full canonical or full-Q222 lift.

## Matrix convention and reproduction

The source order is (d,lambda,n,m,C,A,D,alpha,beta).
The Alice auxiliary order is (tau,p,q). An entry is
[parameter,row,column,rational_coefficient], with zero-based indices.
Mirror every off-diagonal entry once. The source-dependent values used
to prove positivity are witnesses for free auxiliaries; no nonlinear
equations are imposed on the pencil.

With Python 3 and SymPy, run:

    python check_and_extract.py
    python check_native_contact.py

The first command reconstructs both arrays from the square decompositions
and checks all source and auxiliary coefficients against the homogeneous
four-state parent. The second verifies the native contact zero exactly.
The analytic source estimates are proved in the PDF. Compile the TeX
source twice with pdfLaTeX.
