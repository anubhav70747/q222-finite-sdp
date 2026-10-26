REPRODUCING THE EXACT ALGEBRA CHECKS

Python 3 and SymPy 1.14.0 are required. In a suitable Python environment:

  python -m pip install sympy==1.14.0
  python checks/run_checks.py

Run from the proof-project directory. The driver writes fresh outputs and
logs under checks/recomputed and preserves the original checks/receipts.
The programs can take several minutes. No diagnostic program is required
to compile or read the written proof.

The seven main suites contain 340 named exact checks. The last six names
in the weak-C2 suite are the new bounded-a checks; the small separate BA
receipt is a subset, not six additional checks. The three supplemental
support tables each enumerate 128 Boolean event-presence assignments.

The scripts independently transcribe the relevant mathematical formulas.
Portable copies differ from the audited scripts only by removal of an
absolute path to the temporary SymPy installation. The support-table script
compares the printed NC14 conditions with the expanded coefficient receipt;
after a full rerun it uses the freshly recomputed coefficients.

These diagnostics verify selected algebraic identities. They do not prove
inequalities, source coverage, separated tensor membership, or the global
theorem automatically. Those arguments are in the proof and review reports.
