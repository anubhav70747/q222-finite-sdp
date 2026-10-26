# Complementary physical C1 cells: independent bounded audit

Date: 26 October 2026. Source pinned: `8960b3fead3ce6b185bb1387d148fb87d05d3d23`, `latex/certificates/c1_short.tex`. Scope excludes the separately assigned BN/UN, CA/SE, EL/NU, and independently replayed FM cells. No global C1 or Q222 completion is claimed.

## Exact new replay: NC6 and NC7

Put `K=1+x`, `R=1+x*r`, `T=1+x*h`, `G=h+z`, `theta=1-h`, `E=1-x*eta*h`, `y=x*eta*z`,
`p0=1+2*x*h+r*x^2*G+x*z`, `p1=x*(K*z+x*h*theta)`, `Pd=p0+eta*p1`, and
`B5=K*(1+y)*(r+eta)-(eta*T-1+r)*Pd`, `Bb=K*x*(r+eta)*G-theta*Pd`.
Set `Acal=theta*T+K*z`, `Ccal=R*G*(1+r*x*h)` and
`Bcal=h^2*r*x^2+h*r*x*(2-h)+h^2*x+h*r*x^2*z+h*x*z+h*x+h+r*x*z+x*z+2*z+1`.
Then exact polynomial expansion yields

`h*B5=-Acal*E^2+Bcal*E-Ccal`.

For
`a2=K*x*z-T*p1`, `a1=K*(1+x*z*r)-T*p0+(1-r)*p1`, `a0=K*r+(1-r)*p0`,
`b1=K*x*G-theta*p1`, `b0=K*x*r*G-theta*p0`, exact elimination gives

`b1^2*B5-a2*Bb^2-(b1*a1-2*a2*b0)*Bb=b1^2*a0-b1*a1*b0+a2*b0^2`
`=x^2*h^2*K^2*R*G^2*[1+2*x*h+r*x^2*h^2+(1-r)*x*z]`.

These were checked as polynomial equalities (not at selected parameters). On native physical domain `0<x,r,h<=1`, `theta<2*z`, `B5>=0`, the first identity implies `E>=1/18`. If x, y=x*eta*z are positive bounded units, the native gap `x*eta*h<1` and `x*z*r<1` force eta,z to be bounded units. Then `h` is a unit by `Bb|theta=1=-1-x*z*(1-r)` and `Bb>=0`. All terms on the NC7 RHS have positive unit residue. Consequently B5 and Bb cannot both be infinitesimal. The exact NC4 formula then makes at least one of D/d, beta/d strong, while C,A,m,alpha are strong. Thus five consecutive event weights are bounded below by positive reals and the *separately scoped* five-event theorem applies.

## Complementary chart ledger, proof status

1. k unit and eta,zeta not both infinitesimal: the unit/unit contact shortcut above is directly supported; UI and UH remaining shape regimes require the source's CU1–CU13 / CR1–CR13 and UH1–UH20 full-state allocations. These were read but not independently reconstructed line-by-line here.
2. k infinitesimal, eta,zeta infinitesimal: UJ1–UJ5 plus the bounded reciprocal chart UE1–UE4 are source statements for all e; these depend on exact EZ and IW full-parent certificate bodies and their projective joins. Not independently re-proved here.
3. k,zeta infinitesimal, eta unbounded: SC1–SC20 (q-unit), QI1–QI9 and QR1–QR3 (q-infinitesimal) provide the two cells. Their strict Gram conclusions and all-angles joins remain source-derived dependencies, not a consequence of NC7.
4. k,zeta infinitesimal, eta unit: BE1–BE7 for unit r,e and the bounded-reciprocal/centered-A CA cell are separate. CA/SE is outside this assignment. For e unbounded with r unit the source relies on the exact critical-current and contact conductor (UE3); not independently reverified here.
5. k infinitesimal, zeta noninfinitesimal: with L=x*zeta unbounded, unit eta uses the bounded source frame in `pd4b92bf7de71`, unbounded eta uses the separate inverse-Alice theorem; L unit eta-small and eta-unit use LC/LU strong-II adapters, eta-unbounded uses the bounded-z theorem; L small with eta bounded uses SB, IC/JO/RO and simultaneous r,j frame. L small with eta unbounded has five kappa,nu,ell cells: FM and BN/UN are separately assigned; other sectors retain the component references in the pinned source.

The source partition's union is a logical exhaustive partition, but this is not the same as revalidating every imported full-parent SOS identity. In particular none of the results above establishes the noncircular complete complementary C1 all-tests theorem by itself.

## First unresolved implication for this audit

For each residual complementary cell named above, prove that after exact removal of original positive packets, the **complete scalar, both first normals, and hard quadratic block** admit a bounded coefficient Gram with positive definite real residue on every remaining independent projective testing chart. It is not enough to quote the stated strict residue, an affine coefficient identity, or a partition label; the source-specific monic-division coquotients and caps must have opening-uniform bounds and all testing joins must be checked. The first not independently closed interface in the k-unit route is the CU1–CU13 / CR1–CR13 full-matrix joined Gram across original opening limits; the audited NC7 shortcut does not cover UI/UH with one infinitesimal actual contact half-angle.

No GitHub publication was performed; this note is suitable for independent coordinator review before publication.
