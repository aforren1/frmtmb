# Lane sampfix last round: R1, R2, R3 in the findings.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, ":\n", old)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("## Version floor",
"## Last round (re-check R1 to R3)

- **R1, fixed.** `draws_laplace_watch()` spent nearly all its time in
  `unlist(v)` building one name per cell. It now uses
  `unlist(v, use.names = FALSE)`. Measured with the reviewer's
  `dev/sampfix-rev-r2-03-cost.R` (`y ~ x + (1 | g)`, 20000 rows in 200
  groups, 1000 draws, data seed 5, draws seed 1, `re_formula = NA`,
  3 repetitions, medians), laplace over full draws:

  | run | `posterior_predict()` | `posterior_epred()` | watch alone, 1000 calls |
  |---|---|---|---|
  | before (`13-cost-before.txt`) | 1.964 | 2.000 | 3.670 s |
  | after (`13-cost-after.txt`) | 1.037 | 1.050 | 0.160 s |
  | after, again (`13-cost-after2.txt`) | 1.067 | 1.024 | 0.130 s |

  The machine was loaded during the before run (full-draws
  `posterior_predict()` 5.56 s against 3.28 s after), so read the ratios
  and the watch-alone time, not the absolute seconds; the reviewer's
  own before run gave ratios 1.946 and 1.794 and 2.160 s. Results are
  `identical()` between laplace and full draws in every run.
- **R3, fixed.** The watch compared `!is.finite()` patterns, so a later
  draw that overflowed to `Inf` on its own was refused as a read. It
  compares `is.na()` patterns now. Checked first on every watched path
  (`dev/sampfix-12-fillpaths.R`, `dev/sampfix-log/12-fillpaths.txt`,
  data seed 77, draws seed 1): with the probe off and the watch replaced
  by a recorder, every read of the `NA` fill came out `NA` or NaN and
  none came out `Inf`, on 38 combinations. The paths were
  `posterior_epred`, `posterior_linpred` and `posterior_predict`'s
  dpars, in sample and at `newdata`, on gaussian, poisson, bernoulli,
  lognormal with `(1 | g)` in mu and sigma, cumulative, a smooth, and a
  nonlinear `exp(a)^k` body; plus `conditional_effects()` at its default
  and at `re_formula = NULL` on a smooth plus `(1 | g)`, the `mi()` model,
  `pp_mixture()` and `posterior_epred()` on a grouped mixture. The
  `hypothesis()` path had no read to observe (an `sd_` needs none). A
  test (`test-laplace-draws.R`: lognormal, `b_Intercept = 720` on draw 3)
  was seen to fail on the build before the change
  (`dev/sampfix-log/r3-laplace-before.txt`) and passes now.

## Version floor")
rep1("   above is the same code path seen from `conditions`.",
"   above is the same code path seen from `conditions`.
7. **R2 (reviewer), a NaN the watch can still pass.** In
   `bf(yn ~ log(c1) * x + exp(a)^k, c1 ~ 1, a ~ 1 + (1 | g), k ~ 1,
   nl = TRUE)`, with `c1 = -1` AND `k = 0` on draw 1 only, draw 1 is NaN
   on every row for its own reason, so the probe sees NaN at both fills
   and the watch takes the all-`NA` pattern as its reference; draws 2 to
   6, which read `a`'s group effects, match it. `posterior_epred()`
   returns 480 of 480 cells non-finite where the full draws give 80.
   It needs two exact parameter values at one draw and the result is
   NaN, never a finite wrong number. Not fixed; comparing each draw with
   itself at fill 0 would close it at the cost of a second evaluation
   per draw. `dev/sampfix-rev-r2-01-watch.R`, section 1a.")
writeLines(s, f)
cat("ok\n")
