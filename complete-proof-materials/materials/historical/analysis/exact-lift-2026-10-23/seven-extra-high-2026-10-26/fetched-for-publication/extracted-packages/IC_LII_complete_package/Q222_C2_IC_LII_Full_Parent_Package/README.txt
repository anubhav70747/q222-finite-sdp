Q222 C2 weak-IH IC / strong-II LII: full-parent proof package

PROOF.md: statement, exact parent/source/Gram identities and analytic residue admission,
           singular-source endpoint cases and complete projective testing joins.
run_checks.py: sequential exact symbolic replay; requires Python 3 and SymPy.

Checks, all purely symbolic:
  Q222_IC_LII_small_identity_check.py  : IC8/IC8a, LII15(a-e)/LII20, WIE5
  Q222_IC_normal_symbolic_check.py    : both unrestricted weak-IH normal derivatives/monic quotients
  Q222_IC_graph_support_check.py      : the four-monomial IC7 contact residual
  Q222_IH_endpoint_source_check.py    : simultaneous endpoint rho=0,1 source cancellations
  Q222_WIE4_full_graph_identity.py    : removal of full large-D positive square
  Q222_LII_full_parent_readback.py    : LII4/LII6, all four retained packets + full parent
  Q222_LII_coeff_extraction_probe.py  : bounded dependent-root gauge, parities, phi and Pxi support
  Q222_IC_LII_residue_check.py        : IC9, LII19 and remote OA8 determinants

Run: python run_checks.py
There are no numerical sweeps, semidefinite solvers, or approximate checks.
The stated tensor-SOS results concern these two native C2 leaves ONLY; they do not
establish the complete finite SDP or a global uniform feature space.
