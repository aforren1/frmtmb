# Lane seflags: what the curved-ridge correction costs. For each design,
# the fit with the default check against the same fit with
# check_se = "ignore" (the control, the same code in both arms), the
# minimum of `reps` rounds interleaved, plus a count of the objective
# and gradient evaluations the check makes after the fit (counted on a
# wrapped objective, so load-independent). Run once per build.
#   Rscript dev/seflags-cost.R <lane|base> [reps]
a <- commandArgs(trailingOnly = TRUE)
libs <- switch(a[1],
  lane = c("C:/Users/adf44/source/r/wt-seflags-lib",
           "C:/Users/adf44/source/r/rellib-r7"),
  base = "C:/Users/adf44/source/r/rellib-r7")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
reps <- if (length(a) > 1) as.integer(a[2]) else 5L
cat("lib", find.package("frmtmb"), "\n")
designs <- list(
  glmm_uncentered = function() {
    set.seed(1)
    d <- data.frame(x = runif(400, 100, 101), g = factor(rep(1:20, 20)))
    d$y <- rpois(400, exp(-20 + 0.2 * d$x + rnorm(20, 0, 0.4)[d$g]))
    list(f = y ~ x + (1 | g), fam = poisson(), d = d)
  },
  lm_poly5 = function() {
    set.seed(5)
    d <- data.frame(x = runif(200, 1, 2))
    d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
    list(f = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), fam = gaussian(),
         d = d)
  },
  lmm_poly5 = function() {
    set.seed(5)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    list(f = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + (1 | g),
         fam = gaussian(), d = d)
  },
  sleep = function() {
    d <- lme4::sleepstudy
    list(f = Reaction ~ Days + (Days | Subject), fam = gaussian(), d = d)
  }
)
for (nm in names(designs)) {
  s <- designs[[nm]]()
  fitit <- function(ch) suppressMessages(suppressWarnings(
    frm(s$f, family = s$fam, data = s$d,
        control = frmtmb_control(check_se = ch))))
  t <- matrix(NA_real_, reps, 2)
  for (r in seq_len(reps)) {
    t[r, 1] <- system.time(fitit("ignore"))[["elapsed"]]
    t[r, 2] <- system.time(fitit("warning"))[["elapsed"]]
  }
  # the check's own evaluations, on a fit made with the check off
  f <- fitit("ignore")
  n_fn <- 0L; n_gr <- 0L
  fn0 <- f$obj$fn; gr0 <- f$obj$gr
  f$obj$fn <- function(...) { n_fn <<- n_fn + 1L; fn0(...) }
  f$obj$gr <- function(...) { n_gr <<- n_gr + 1L; gr0(...) }
  invisible(suppressWarnings(suppressMessages(ns$se_check(f,
    frmtmb_control(check_se = "warning")))))
  cat(sprintf("%-16s ignore %.3f s | check %.3f s | ratio %.3f | check evals fn %d gr %d\n",
              nm, min(t[, 1]), min(t[, 2]), min(t[, 2]) / min(t[, 1]),
              n_fn, n_gr))
}
