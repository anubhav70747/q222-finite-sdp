# Apply the additive Q222 review overlay

Base: `1ad5125dbfc9c56881c283724b04994e3d50ec00` in `anubhav70747/q222-finite-sdp`.

Inspect these additions in a checkout at the base revision and run:

```text
node extracted-lifts/global-weighted/angular-order13/verify_angular13.js
node independent-review/checks/contact-ray-exposure/verify_contact_seeds.js
node independent-review/checks/l12-half-power/check_dictionary.js
```

Read `reconciliation/2026-10-26/README.md`, the 107-row crosswalk and `BLOCKERS.md`. Existing sources/PDFs are retained. The proposed documentation patch is a separate reviewable file. Preparation does not update any branch.
