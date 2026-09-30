# Lane sampfix nits round: corrections to dev/sampfix-findings.md (S4).
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, ":\n", old)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("Laplace-conditional draw is right to within Monte Carlo error. It was
NOT shipped, for two reasons that need the user: `posterior_epred()`
would depend on the random seed (brms's does only for new levels), and
`log_lik()`, `ranef()`, `coef()` and `loo()` refuse laplace draws by
design, so filling `b` in for the predictive methods alone would make
the methods disagree about what a draw contains. Shipping it means
changing all of them together, with a way to keep the filled `b` fixed
across calls.",
"Laplace-conditional draw is right to within Monte Carlo error on this
gaussian model, where the Laplace approximation of `b` given `theta` is
exact (for another family it is itself an approximation). It is NOT
shipped. **The user decided this on 2026-09-29:** `laplace = TRUE`
samples an approximate posterior whose use is checking the
approximation, and anyone who wants predictions samples without it. Two
further reasons: `posterior_epred()` would depend on the random seed,
and `log_lik()`, `ranef()`, `coef()` and `loo()` refuse laplace draws
by design, so filling `b` in for the predictive methods alone would make
the methods disagree about what a draw contains.")
rep1("just handed. On the reference, 28 calls of the 75 in the sweep died on
`@` for a `stanfit = NULL` object; on the lane, 0.",
"just handed. On the reference, 26 calls of the 75 in the sweep died on
`@` for a `stanfit = NULL` object, and 2 more (`nuts_params()`,
`log_posterior()`) on bayesplot's \"no applicable method\"; on the lane,
0 died on `@`. (The first version of this file said 28, counting those
2 as `@`; corrected in the nits round, `dev/sampfix-10-summary.R`.)")
rep1("  values do not have. So does frmtmb.sample on frmtmb 0.65.0.",
"  values do not have. frmtmb.sample's floor rises to the frmtmb that
  carries this at consolidation (other lanes need new core exports), so
  the names do not depend on which core is installed.")
rep1("reference build at the same seed, the stored draws mapped the lane's
way, `posterior_epred()` and `fixef()` differ by exactly 0, and the
reference's own draws object (old names) read by the lane build gives
the same `posterior_epred()` and `fixef()`, difference 0.",
"reference build at the same seed, the stored draws mapped the lane's
way and `fixef()` differ by exactly 0, and `posterior_epred()` by 0 on
these four. The reviewer's 11 fits (`dev/sampfix-rev-04-ordinal.R`)
found `posterior_epred()` 1 ulp (2.2e-16) apart on 4 of them, all with
cumulative thresholds, where the inverse map `c(tau1, log(diff(tau)))`
does not return the sampled log increments bit for bit; so the claim is
\"within 1 ulp\" for `posterior_epred()`, and exactly 0 for the stored
draws and `fixef()`. The reference's own draws object (old names) read
by the lane build gives the same `posterior_epred()` and `fixef()`,
difference 0.")
rep1("sampler seed 3). 45 of 150 results changed.",
"sampler seed 3). 52 of 150 results changed (45 before the nits round,
which reworded the refusals).")
writeLines(s, f)
cat("ok\n")
