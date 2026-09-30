# Lane sampfix nits round: the nits section, M2/M7 and the floor.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, ":\n", old)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("cumulative thresholds, where the inverse map `c(tau1, log(diff(tau)))`
does not return the sampled log increments bit for bit; so the claim is",
"cumulative thresholds, which is consistent with the round trip through
`c(tau1, log(diff(tau)))` not returning the sampled log increments bit
for bit (not measured separately); so the claim is")

rep1("   draws for `sigma` and now also for the thresholds. Not changed.",
"   draws for `sigma` and now also for the thresholds. Not changed.
5. **M2 (reviewer): a harmless false alarm.** `(0 + x | g)` at `newdata`
   with `x = 0` is refused on laplace draws, because the probe's `NA`
   times 0 is `NA`, although the prediction there does not depend on
   `b`. The user is told to sample without `laplace = TRUE`, which
   works. `dev/sampfix-rev-01-probe.R`, rerun on the final build in
   `dev/sampfix-log/11-probe-lane.txt`.
6. **M7 (reviewer), pre-existing:** `conditional_effects(re_formula =
   NULL)` on draws draws a NEW level from `N(0, Sigma(theta))` (core
   `ce_draw_new_levels()`), where brms 2.23.0 calls its predictions with
   `allow_new_levels = TRUE` and the default `sample_new_levels =
   \"uncertainty\"`, which draws from the existing levels' effects. Full
   draws as well as laplace ones (`dev/sampfix-rev-08-brmsce.R`). Item 1
   above is the same code path seen from `conditions`.")

rep1("frmtmb.sample needs no new core EXPORT. It reads a new family field,
`post$ord_thresholds_raw`, through the existing `brms_fixef_rows()`,
and falls back to the old names without it. Its DESCRIPTION floor can
stay at `frmtmb (>= 0.65.0)`; the ordinal renaming only takes effect on
the core that carries this lane. Both packages: a patch bump would
understate the draws-name change; minor for frmtmb.sample (a visible
rename of draw columns), patch for frmtmb.",
"frmtmb.sample needs no new core EXPORT. It reads a new family field,
`post$ord_thresholds_raw`, through the existing `brms_fixef_rows()`.
Settled at consolidation (M5): frmtmb.sample's floor rises to the next
frmtmb anyway, because other lanes need new core exports, so the draws'
names never depend on which core is installed; NEWS says so. The
fallback to the internal names stays in the code for a family that
declares no inverse. Bumps: minor for frmtmb.sample (a visible rename
of draw columns), patch for frmtmb.

## Nits round (after review, dev/reviews/2026-09-29-sampfix.md)

- **S1.** A laplace refusal suggests `re_formula = NA` only when the
  function the user called takes it and the same call at
  `re_formula = NA` passes the probe on this model
  (`draws_laplace_refuse()`); otherwise it says only \"Sample without
  laplace = TRUE\". Rerun of the reviewer's message script on the final
  build (`dev/sampfix-log/11-messages-lane.txt`): `bayes_R2()`, the
  smooth-only model and `pp_mixture()` carry no hint now.
- **S2.** Every refusal names the function the user called:
  `draws_as_caller()` records the outermost caller, and `fitted()`,
  `predict()`, `residuals()`, `predictive_error()`,
  `predictive_interval()`, `bayes_R2()`, `pp_check()`, `loo()`,
  `waic()`, `psis()`, `loo_compare()`, `hypothesis(scope =)` and
  `coef()` set it. In the rerun every refusal names the called function;
  the reviewer's check flags the second `pp_check()` line only because
  that message begins `pp_check(type = 'loo_pit_overlay')`, which names
  it with its type.
- **S3.** `?frm_sample` leads with the user's reason and no longer cites
  a dev/ script.
- **M1.** `draws_laplace_watch()`: on laplace draws, a draw whose result
  has non-finite cells where the first draw's has none is refused with
  the same message. Compared with \"any non-finite output\", this keeps a
  row that is `NA` at every draw (a missing covariate) from being taken
  for a read. The reviewer's construction (`c0 + exp(a)^k`, `k = 0` on
  draw 1) is a test in `test-laplace-draws.R`, seen to fail on the build
  before this round (`dev/sampfix-log/nits-laplace-prenits.txt`, the
  failure at line 233) and passing now; the probe rerun no longer lists
  it among its FAILs.
- **r_eff.** `test-ppcheck-loo.R` pins brms's rule: the chain structure
  when every draw is used, one chain for a subset, with the two rules
  shown to differ on the same draws.
- **S4.** The corrections above (26, not 28; within 1 ulp);
  `dev/sampfix-log/05-lane.txt` regenerated on the final build (all
  four types, all draws and the subset: max |lw diff| 0, max |k diff| 0).")
writeLines(s, f)
cat("ok\n")
