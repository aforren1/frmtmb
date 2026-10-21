## The 0.73.0 round, 2026-10-20 to 2026-10-21

`dev/round-20261021.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.26
and 0.3.32 (BLAS and LAPACK, 4 threads). Each lane has
`dev/<lane>-findings.md` and `dev/reviews/2026-10-20-<lane>.md`. Four
lanes off `4854e425` (0.72.0), merged in this order:

- `loodraws`: brms's `loo_R2()`, `loo_epred()`, `loo_linpred()`,
  `loo_predict()`, `loo_predictive_interval()` and
  `loo_model_weights()` on draws; brms's `r_eff` rule on a subset of
  draws; `posterior_average()` weights from `loo_model_weights()`; ten
  generics in core (`R/loo-generics.R`) with their refusals on a fit.
- `cdfs`: the CDFs of 23 more families (`R/families-cdf.R`), the
  log-scale window code of `row_lpdf()` and the `cdf_log` slot, brms's
  reading of a discrete censoring bound, brms's two `trunc()` checks,
  the exgaussian closed form, `simulate()`'s inverse-CDF fallback.
- `nb2dm`: `brmsfamily("negbinomial2")` and `dirichlet_multinomial()`,
  and the log-scale Dirichlet draw.
- `fuzzfix`: quadrature refuses a row in more than one block (the
  exact iterated route is `dev/fuzzfix-nested-quad.patch`), and the
  refusal of a weighted likelihood with no maximum
  (`check_weights_unbounded()`).

Versions: frmtmb **0.73.0** (minor, breaking); frmtmb.sample
**0.21.0** (minor), floor frmtmb 0.73.0 (measured: 0.21.0 does not
install over 0.72.0, `dev/rel073-floor.sh`). The six other extensions
keep their versions and floors: no lane changed a file of theirs. The
release library is `C:/Users/adf44/source/r/rellib-r11`; `rellib-r10`
keeps 0.72.0. In brms's ported suite bins 1 and 2 pass 607 of 823
(73.8%; 574 at 0.72.0), 607 of 706 (86.0%) without the decisions; bin
1 alone 392 of 494, unchanged.

### Decisions made at the consolidation, for the user to confirm

- **The emmeans leak is fixed in the generator**, with
  `withr::local_package("emmeans")` at the file's top level, which
  detaches at the end of the file (measured; `teardown_env()` would
  not).
- **rstan's Rcpp module classes are tagged as a load** in
  `dev/release/run-session.R` ("[with load of a Stan model (rstan)]"),
  measured to come from rstan alone; every other global stays
  untagged.
- **`loo_predict` joins core's `_pkgdown.yml` index**, which lane
  loodraws had left out (`pkgdown::check_pkgdown()` failed on it).
- **The NEWS sentences over 25 words are split** (28), as at 0.72.0.
- **`?posterior_average` links its weights to `sample-loo_model_weights`**
  by name: roxygen had resolved the bare `[loo_model_weights()]` to
  `loo::` or to the local re-export, depending on the library it ran
  against.

The release review (`dev/reviews/2026-10-21-release.md`, SHIP AFTER
FIXES) is applied ("After the release review" in the record): the
weights and asym_laplace check now runs on the uncensored rows of a
`cens()` model and carries `zero_inflated_asym_laplace()`; the four
core redirects of `kfold`, `reloo`, `loo_moment_match` and
`loo_subsample` point at `reference/loo_predict.html`; lane adhess's
filed poisson `trunc(lb = 1)` seed-19 item is closed by lane cdfs's
log-survivor normalizer (measured by the review).

### What the 0.73.0 round is evidence for

**A lane that adds an Rd file to another package must add its site
index entry.** Lane loodraws moved ten generics into core with a new
`man/loo_predict.Rd`; its core suite and R CMD check passed, and only
`pkgdown::check_pkgdown()` (the site workflow's own check) stopped on
the missing topic. A lane that adds or renames a topic runs
`check_pkgdown()` for the package it touched.

**A `teardown_env()` cleanup is a run's, not a file's.** The obvious
fix for an attached package, `withr::local_package(..., .local_envir =
teardown_env())`, keeps it attached for every later file of a session;
only the file's own environment detaches it at the file's end
(`dev/rel073-local-package.R`). The one-session gate is the instrument
that shows the difference.

**A three-way numbering clash is cheap only if one pass renumbers.**
Three lanes each appended a "brms-30". Renumbering one map at a time
(30 to 31, then 31 to 32) would have shifted cdfs's entries twice;
`dev/rel073-resolve.R` maps each text once, and two passages of a
review that discuss the numbering itself keep their numbers with a
note.
