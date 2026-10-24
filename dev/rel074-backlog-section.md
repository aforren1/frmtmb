## Filed at the 0.74.0 release (2026-10-23)

Three lanes off 0.73.1, the first round under the core boundary rule
(`dev/round-20261023.md`, `dev/core-boundary.md`). What each lane, its
reviewer or the consolidation filed and did not fix, with the record to
read. Nothing was posted upstream.

### frmtmb.survival failed R CMD check --as-cran: FIXED at the release

- **`test-saved-fits.R` failed under R CMD check** (7 failures, the
  whole "not installed" arm; `Status: 1 ERROR, 1 NOTE`;
  `dev/rel074-check/survival-before/`). The test built a library stack
  without frmtmb.survival by dropping every library that holds it.
  Under R CMD check --as-cran the tests see one temporary library that
  holds the package and all its dependencies, then R's own
  (`dev/rel074-log/survival-check-probe.txt`), so the child process had
  no `generics` and could not load frmtmb. FIXED at the release, as a
  release fix the coordinator approved (`dev/rel074-release-fix.R`):
  both packages copy one helper, `helper-absent-package.R`, with the
  two stacks and frmtmb.spatial's probe; frmtmb.survival gained the
  "installed but unloadable" arm, which runs under R CMD check, and
  only the "not installed" arm skips there. Check after the fix:
  `Status: 1 NOTE` (V8), FAIL 0, SKIP 4, PASS 366
  (`dev/round-20261023.md`, "R CMD check --as-cran").
- Still open from it: a lane that writes a test that builds its own
  library stack runs R CMD check on the final test. The two copies of
  `helper-absent-package.R` have no check that they stay identical, as
  `helper-scale.R`'s copies have in `scale-helper-identical.yaml`.

### Stale text found at the consolidation

- `.github/workflows/check-frmtmb-learn.yaml` compares itself with "the
  spline job" twice (comments); frmtmb.spline's job is gone.

### brms spellings of `cox()`

- **brms's `bhaz()` addition term.** `time | bhaz(df = 6) ~ x` (with
  `knots =` and `Boundary.knots =`) is refused by name since 0.74.0,
  with a pointer to `cox(df =, degree =, intercept =)` in
  frmtmb.survival, but it is not accepted. frmtmb.survival could take
  `bhaz(df =, knots =, Boundary.knots =)` into `cox()` and refuse
  `bhaz(gr = )` by name, since frmtmb has no stratified baseline.
  Record: `dev/splitsurv-findings.md` ("Defects found, not fixed");
  `dev/reviews/2026-10-23-splitsurv.md`, finding 7.
- **`brms::cox(bhaz = list(df = 6))` passed to `frm()`** stops with
  "unused argument (bhaz = list(df = 6))", which names neither the
  frmtmb spelling (`cox(df = 6)`) nor a reason. Same records.

### Coefficient names

- **The coefficient of `exp(x)` is named `.frm_expx`; brms names it
  `expx`.** `fixef()` gives `Intercept, .frm_expx` and `variables()`
  `b_.frm_expx`, on 0.73.0 (rellib-r11) and on 0.74.0 alike. brms
  2.23.0's `standata()` has the columns `Intercept expx` and
  `default_prior()` the coef `expx`. The parser's placeholder for a
  special name (`sub_specials()`, `dev/rtmb-pitfalls.md` item 9) leaks
  into the brms names, probably for every function it rewrites this
  way. It breaks `set_prior(coef = "expx")` and `hypothesis()` on
  brms's spelling. Found by the spatial review
  (`dev/reviews/2026-10-23-spatial.md`, "A pre-existing name defect,
  outside this lane"); not caused by this round.

### The covariance-structure seam

Record: `dev/spatial-findings.md` ("What I did not do") and
`dev/reviews/2026-10-23-spatial.md`.

- frmtmb.sample's `draws_cov_fast()` keeps its own list of field
  structures. A registered field is correct there (the map falls back
  when it disagrees with `hyp_env_vals()`) but slower. A core reader
  export would close it.
- `rr()`'s factor-space expansion, `gp()` and `hsgp()` as specials, and
  the Student-t sharing of `us_t` and `diag_t` are not on the published
  entry shape.
- The weights check (`unbounded_full_rank_blocks`) fails open for a
  registered structure, and autoscale skips one without a message
  (the spatial review's table of consumers). Neither is wrong for the
  built-in structures.
- The refusal of `sar()` and `fcor()` says "The supported residual
  correlation terms are ar(), ma(), arma(), cosy(), unstr() and car()"
  (`R/parse.R`). `car()` is in frmtmb.spatial since 0.74.0, so the
  sentence should say so. Found at the consolidation.

### A saved fit that needs two extensions

- A fit with a family from one extension and a structure from another
  (`cox()` with a `car()` term) is refused at the first method that
  needs a missing package, by that package's name, and never with
  "could not find function" (`dev/rel074-saved-mixed.R`,
  `dev/rel074-log/saved-mixed.txt`). But the refusal names one missing
  package at a time: without both, `summary()` names frmtmb.survival
  and `print()` names frmtmb.spatial, and the user learns of the second
  only after installing the first. One message could name both. Found
  at the consolidation.

### Releases in lockstep

- frmtmb.latent 0.6.2 does not load with frmtmb 0.74.0: its `.onLoad()`
  registers `hmm()` rules that name the `mixture_mvn` feature, which
  core no longer has. The release is in lockstep (frmtmb has no
  existing users). From now on an extension can declare a feature that
  another package supplies with `expects =`, so a later move does not
  break an older extension's load. Record: `dev/mvnlatent-findings.md`
  section 6; `dev/reviews/2026-10-23-mvnlatent.md`, finding 5.
- frmtmb.sample suggests frmtmb.latent with no version. Its
  `test-simulators.R` block "mixture_mvn simulates its class
  covariances" fails, not skips, with a frmtmb.latent older than 0.7.0.

### mixture_mvn() and the latent tests

Record: `dev/mvnlatent-findings.md` section 7 and fix round 1;
`dev/reviews/2026-10-23-mvnlatent.md`, finding 11.

- The moved tests keep their absolute tolerances (`1e-5`, `1e-4`,
  `1e-10`, `0.3`), so that the move stayed verbatim. They should be
  ratios to a quantity the run measures (`dev/lane-rules.md`).
- frmtmb.latent's `test-prior-dpar.R` reaches three core internals
  (`frmtmb:::prior_design`, `resolve_priorlist`, `par_template_names`)
  through its `pd_reached()`, and the moved block depends on it too.
- `predict(newdata =)` refuses a `mixture_mvn()` fit with "its draws are
  structured (a hidden state sequence, a group-level latent class, or a
  correlated residual)" (`R/predict-brms.R`): none of the three is a
  multivariate normal mixture. Pre-existing.

### Smaller items from the splitsurv review

Record: `dev/reviews/2026-10-23-splitsurv.md`, finding 8.

- `set_prior("dirichlet(2)", class = "Intercept", dpar = "sigma")` is
  refused when it is written, while the same prior written
  `class = "sigma"` is refused at the fit. brms refuses neither.
- `print()` of a curve refuses `digits =` with the generic "has no
  argument" message. Users of data frames will try it; the message
  could name it.
