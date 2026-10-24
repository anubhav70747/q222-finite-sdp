Puiseux-Gram lift reproduction. Run: python gram_lift.py --k K (--e "lam,n,m,C,A,D,alpha,beta" | --case NAME) [--robust]
Cases (cases.json): e_sep, sym36, C1_7, plus the 9 k=2 gaps random8/22/77, C1_1/11/14, C2_2, cbar, cdag. test_output.txt: sample runs.
gaps_k2.json: the 9 confirmed k=2 gaps (hp.py brackets): fam, idx, e (lam,n,m,C,A,D,alpha,beta), ds (=d*), ratio, status, m_lo/m_hi (k=2 margins at 0.999/1.001 d*), final (hp ratio).
sweep6.jsonl (one case per line, k=6, 302 cases): fam, idx, e, ds (d*), ratio (Clarabel d_6/d*; noisy, can be >1), status, sizes (Gram blocks),
  t_lo/ts_lo/pres_lo (margin t at 0.999 d* in the equilibrated basis, status, equality residual), flagged, B (delta-margin results for flagged cases: t, tstat, pres per delta), time.
gap6r_k6.jsonl: k=6 delta-margins on the 9 gaps: t[delta], tstat[delta], pres[delta] for delta in 1e-2,1e-3,1e-4; ref_t centering margins; bscale_range.
Reading t: t>0 and t/delta constant across delta => feasible down to (1-delta) d*; t<0 => gap. t is relative to the equilibrated Gram diagonal, so its size is only comparable within one method.
Clarabel max-d ratios are unreliable on face (zero event) and multiscale cases; use the delta-margin. symcheck.py: sign-symmetry check (prints 0, 0).
Deps: numpy, scipy, sympy, mpmath, clarabel. Reference values: e_sep k=2 ratio 1.0000015, k=6 1.0000036; random8 k=2 gap (hp 0.79035; Clarabel 0.80-0.81), k=6 t>0.
