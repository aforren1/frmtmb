### Filed by lane illse (`dev/illse-findings.md`, sections 9 to 11)

#### Open - medium

- **Degree 7, seed 4 of the 80 sd-at-zero fits now loses its
  coefficient SEs** (`dev/illse-log/fx/r2-edge4-lane.txt`, `-ob32.txt`).
  A larger group sd gains 1.1e-3 there, above `grad_tol`, so the fit
  warns that it may not have converged and keeps `sdreport()`'s
  standard errors, which are NaN. 0.70.0 named theta_1 "short" and kept
  finite coefficient SEs. Safe (a NaN with a warning, not a wrong
  number), but theta_1's verdict now equals 0.70.0's on 77 of the 80
  fits, not 78. The coordinator filed it at the consolidation.
- **The influence measures and the D1 and D2 tests** of `anova()` on
  `frm_multiple()` fits (`cooks.distance()`, `solve(V)`, `solve(U)`)
  read the whole matrix; at degrees 7 and 8 they stop with an error.
- **marginaleffects and a caller's `a' vcov() a`** cannot be made
  accurate on these designs in double precision; they get NaN or the
  `frmtmb_cov_products` warning.

#### Open - low

- Fits with more than 100 outer parameters, a bound that holds a
  parameter or an sd at its edge keep 0.70.0's verdicts (at degree 8
  with `(1 | g)` and the sd at zero only theta_1 is lost on 16 of 20
  seeds; on the others the route does not apply).
- REML, and nonlinear formulas, keep the "illcond" verdict.
- frmtmb.sample's jitter is scaled by the diagonal SEs only, not by
  the factor.
- The cost: 1.45 to 1.66 times 0.70.0's on n = 1e5 with 100 outer
  parameters, and 1.47 to 1.52 with `(1 | g)` (section 10, item 7).
