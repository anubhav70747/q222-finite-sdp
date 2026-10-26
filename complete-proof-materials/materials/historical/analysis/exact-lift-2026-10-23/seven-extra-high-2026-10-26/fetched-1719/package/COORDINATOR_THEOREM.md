# C1 transport in one fixed sixth-root-multiplier cone

**Conditional theorem.** If all strict oriented physical C1 sources `0<k<=1` have PSD `Z16` Grams and all strict physical C2 sources have PSD `Z6` Grams, then the **same one** fixed, globally sound, closed, rational-data SDP cone

`K = {c>=0 : exists Q>=0, (1+v^2)SG2(c) = Z_{6,12}^T Q Z_{6,12}}`,

with `Z_{6,12} = {u^(i/6)v^(j/6)psi_l:0<=i<=6,0<=j<=12,l=W,X,Y,Z}` (364 rows), admits every strict C1/C2 contact ray, including C1 sources initially having `k>1`.

**Exact reflection.** `sigma(d,lambda,n,m,C,A,D,alpha,beta) = (d,n,lambda,C,m,alpha,beta,A,D)`; `G(v)psi=(vW-X,-W-vX,Z-vY,Y+vZ)`; `G^2=(1+v^2)I`, and

`u^2 F_{sigma c}(1/u,v;G(v)psi)=(1+v^2)F_c(u,v;psi)`.

For strict C1 `(k,f,g,a,t)`, set `kappa=k^2`, `eta=f^2`, `zeta=g^2`, `r=1-a eta`, `h=1-t zeta`, `H=r+eta`, `G0=h+zeta`, `M=1+kappa r-t`, `J=1+kappa h-a`. Then

`k'=1/k, f'^2=kappa zeta H/G0, g'^2=kappa eta G0/H, a'=M/(kappa H), t'=J/(kappa G0)`.

All gaps remain strictly positive. If original native `C=2`, reflected physical C weight is `m`, so its oriented native representative is `(2/m)sigma(c)`; homogeneity returns `sigma(c)` without a nonlinear source-cone map.

**Explicit fixed row transport.** `h16=u Z16(1/u,v;G(v)psi)` is a fixed 16-column whose monomial expansions lie in `Z_{6,12}` (exact 32-entry signed sparse map in the JSON). If `F_{sigma c}=Z16^T Q Z16`, then `Dv F_c=h16^T Q h16`. If `F_c=Z6^T Q Z6`, then `Dv F_c=Z6^T Q Z6+(v Z6)^T Q(v Z6)`. Thus both assumptions map into one cone K.

**Closure, dual and assembly.** K is closed by a 364-row positive feature-moment trace bound; `Dv>0` ensures soundness and projective limits. Published singular-density then includes all third-omission C1/C2 contact rays. SG12/SG13 face cones and event rays join with finitely many physical symmetry copies, using the previously proven finite closed-sum theorem. In coefficient notation `T c=coeff(DvF_c)`, `A Q=coeff(Zext^TQZext)`, the exact source dual requires

`b(q)=T^T y+r, A^*y>=0, r>=0`, where `b(q)=(H(q),4p_lambda(q),...,4p_beta(q))`.

The full `Z6` strictly feasible rational source extends to a strictly positive extended Gram because the `Z6` and `vZ6` subgrids together cover all 364 rows, proving attained duality. Neither `A^*y>=0` alone nor unchanged-`Z6` membership suffices.

**Still unproved:** oriented strict C1 Z16 PSD and strict C2 Z6 PSD. The published unoriented C1 Z16 counterexample does not affect the conditional theorem. The previously reviewed qualitative existence proof is not reopened. No Git write occurred.
