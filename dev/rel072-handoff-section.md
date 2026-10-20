## The 0.72.0 round, 2026-10-14 to 2026-10-20

`dev/round-20261020.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.26
and 0.3.32 (BLAS and LAPACK, 4 threads). Each lane has
`dev/<lane>-findings.md` and `dev/reviews/*-<lane>.md`. Nine lanes off
`47cf659d` (0.71.0 plus the test-order fix), merged in this order:

- `testiso`: test isolation (`local_frmtmb_registries()`, the
  compatibility contract with all seven extensions loaded), the
  one-session release gate (`dev/release/run-sessions.sh`), and the
  coverage workflow (core alone, plus one job per extension). No
  package code.
- `catmulti`: brms's `categorical()` and `multinomial()` (`refcat`,
  `refcat = NA` with the reference model's degrees of freedom, the bare
  `multinomial()`), `fitted()` and `posterior_linpred()` of every
  location predictor, `logistic_normal()` in `mvbf()`, the cumulative
  tied-threshold `NaN`.
- `prilabel`: per-coefficient `sd` priors, brms's labels and order of
  nested and multi-membership groups (and of the theta parameters), one
  Student-t latent scale per factor.
- `sampleprior`: brms's default priors in frmtmb.sample, `seed =`,
  `sd_`/`cor_` draw columns; core `uniform()`, `lognormal()`, support
  guards on the tape, the logit Jacobian.
- `bin2port`: bin 2 of brms's suite on 0.71.0's families, with value
  checks (dev only).
- `ndaterm`: every addition term on `newdata`, coded as the fit coded
  it; mixture and `hmm()` term rows; the double softmax of
  `predict(newdata =)`.
- `illse2`: Cook's distance, D1 and D2 through the covariance factor;
  REML and nonlinear bodies in the design's coordinates; bounded
  Newton steps; the 400-parameter cap.
- `adhess`: the exact AD Hessian for every fit without random effects;
  RMS-scaled steps with them; "nocurv"; frmtmb.coupling's refusal.
- `drop1raw`: `drop1()`, `step()`, `update()` and `influence()` refit on
  the data the fit holds, with the fit's bases.

Versions: frmtmb **0.72.0** (minor, breaking); frmtmb.sample **0.20.0**
(minor), floor frmtmb 0.72.0; frmtmb.eam **0.13.0** (minor), floor
frmtmb 0.72.0 and Suggests frmtmb.sample (>= 0.20.0); frmtmb.latent
**0.6.2** (patch), floor frmtmb 0.72.0; frmtmb.coupling **0.6.1**
(patch), floor kept at frmtmb 0.61.0 (measured: its suite passes on
0.71.0; update the two together). Each floor was measured by
installing the extension alone over rellib-r9 (`dev/rel072-floor.sh`).
frmtmb.spline, frmtmb.learn and frmtmb.ode keep their versions. The
release library is `C:/Users/adf44/source/r/rellib-r10`; `rellib-r9`
keeps 0.71.0. In brms's ported suite bins 1 and 2 pass 574 of 823
(69.7%; 549 at 0.71.0), 574 of 706 (81.3%) without the decisions; bin
1 alone 392 of 494 (391).

### Decisions made at the consolidation, for the user to confirm

- **frmtmb.coupling keeps its floor** (frmtmb >= 0.61.0) where the
  brief and lane adhess named 0.72.0: installed alone over 0.71.0, its
  11 test files pass. A floor cannot stop coupling 0.6.0 meeting core
  0.72.0, so its NEWS lead says to update the two together.
- **frmtmb.eam suggests frmtmb.sample (>= 0.20.0)**: lane sampleprior's
  `test-sample-ndt-default.R` in frmtmb.eam fails with frmtmb.sample
  0.19.0's defaults.
- **The two consolidation fixes** (`dev/round-20261020.md`, "Found at
  the consolidation"): `confint_varcorr()` no longer spreads a
  "nocurv" NaN to rows that do not move along it, and lane
  sampleprior's draw-name test takes prilabel's `variables()` order.

### What the 0.72.0 round is evidence for

**A lane's reruns after a fix round are a subset, and the subset can
miss its own regression.** Lane adhess's full suites passed before its
fix round 1; the fix (`sdr_nan_nonpos()`, a NaN row for a non-positive
variance) was rerun on the files the lane chose, and
`test-review-v29.R` was not among them. On the merged build it failed
10 expectations, and on lane adhess's build alone too
(`dev/round-20261020.md`, "Found at the consolidation"). A fix round
that changes package code reruns the whole suite, as the round's
first build did.

**A test that writes another lane's order down is a merge conflict
git cannot see.** Lane sampleprior's new test pinned `variables()` in
0.71.0's order as "brms's names"; lane prilabel changed that order to
brms's. Each lane's suite passed; the merged one did not. When a test
asserts an order or a label, run it against the brms call it claims
to mirror rather than against the build of the day.

**Measure a floor; do not take it from the lane.** Lane adhess said
frmtmb.coupling needs frmtmb 0.72.0; installed alone over 0.71.0, its
whole suite passes, so its floor stays (`dev/rel072-floor.sh`). The
other three extensions whose floors rose fail on 0.71.0, as their
lanes said.
