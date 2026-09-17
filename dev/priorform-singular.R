# Is the fit at test-prior-compat.R:505 genuinely singular on its data?
# The data generator gives g a random intercept (sd 0.7) and nothing
# else: no x slope varying by g, and no h effect at all. So the model's
# (x | g) slope sd and the whole (z | h) block have true value zero.
#   Rscript dev/priorform-singular.R ref|lane
args <- commandArgs(trailingOnly = TRUE)
lib <- switch(if (length(args)) args[1] else "ref",
  ref = "C:/Users/adf44/source/r/rellib-r3",
  lane = "C:/Users/adf44/source/r/priorform-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

make <- function(seed, with_h) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:20, 15)), h = factor(rep(1:15, 20)))
  mu <- 1 + 0.5 * d$x - 0.3 * d$z + stats::rnorm(20, 0, 0.7)[d$g]
  if (with_h) {
    # the two blocks the model asks for, each with a real intercept and
    # slope sd, drawn AFTER the original stream so the seed-77 x, z and g
    # intercepts are unchanged
    mu <- mu + stats::rnorm(20, 0, 0.4)[d$g] * d$x +
      stats::rnorm(15, 0, 0.6)[d$h] + stats::rnorm(15, 0, 0.4)[d$h] * d$z
  }
  d$y <- stats::rnorm(n, mu, 1)
  d
}

run <- function(d, label) {
  w <- character(0)
  fit <- withCallingHandlers(
    frm(bf(y ~ x + z + (x | g) + (z | h)) + gaussian(), data = d),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  vc <- VarCorr(fit)
  sds <- unlist(lapply(unclass(vc), function(b) sqrt(diag(b))))
  cat(sprintf("%-28s warnings=%d  %s\n", label, length(w),
              paste(w, collapse = " | ")))
  cors <- vapply(unclass(vc), function(b) stats::cov2cor(b)[1, 2], 0)
  cat("   sd:", paste(sprintf("%s=%.3g", names(sds), sds), collapse = " "),
      "\n   cor:", paste(sprintf("%s=%.6f", names(cors), cors),
                         collapse = " "), "\n")
  # lme4 on the same data is an independent reading of whether the ML
  # point sits on the boundary of the covariance space
  if (requireNamespace("lme4", quietly = TRUE)) {
    m <- suppressWarnings(suppressMessages(
      lme4::lmer(y ~ x + z + (x | g) + (z | h), data = d, REML = FALSE)))
    cat("   lme4 isSingular:", lme4::isSingular(m), " cor:",
        paste(sprintf("%.6f", vapply(lme4::VarCorr(m), function(b)
          attr(b, "correlation")[1, 2], 0)), collapse = " "),
        " logLik frmtmb - lme4:",
        sprintf("%.3g", as.numeric(logLik(fit)) -
                  as.numeric(stats::logLik(m))), "\n")
  }
  length(w)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("## the test's own data, seed 77\n")
run(make(77, FALSE), "seed 77, as in the test")
cat("\n## seed 77 under the generator the test now uses\n")
run(make(77, TRUE), "seed 77, new generator")
cat("\n## same generator, 10 seeds: how often does it warn?\n")
nw_old <- vapply(1:10, function(s) {
  run(make(s, FALSE), paste("old seed", s)) > 0
}, TRUE)
cat("\n## generator with the effects the model asks for, same 10 seeds\n")
nw_new <- vapply(1:10, function(s) {
  run(make(s, TRUE), paste("new seed", s)) > 0
}, TRUE)
cat("\n<!-- priorform-singular:begin -->\n")
cat(sprintf("old generator: %d of 10 seeds warn; new generator: %d of 10\n",
            sum(nw_old), sum(nw_new)))
cat("<!-- priorform-singular:end -->\n")
