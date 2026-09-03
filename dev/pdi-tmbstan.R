# dev/pdi-tmbstan.R -- does tmbstan sample the tape it was handed?
#
# dev/pdi-repro.R proved the augmented objective is correct and identical on
# both platforms, so the drop happens between MakeADFun() and the draws. This
# script removes frmtmb: two RTMB objectives that differ ONLY by a tight
# normal prior term, sampled with the same seed. The two posteriors cannot
# overlap unless the prior was ignored.
#
# tmbstan() replaces the R closures with the raw tape pointer:
#     mod@ptr <- obj$env$ADFun$ptr ; mod@DLL <- obj$env$DLL
# unless debug = TRUE, which keeps the closures. Running both tells us which
# side of that swap loses the prior.
#
#   Rscript dev/pdi-tmbstan.R

options(warn = 1, digits = 12)
suppressPackageStartupMessages({
  library(RTMB)
  library(tmbstan)
})

cat("##### PDI tmbstan probe\n")
cat("R          : ", R.version.string, " ", R.version$platform, "\n", sep = "")
for (p in c("RTMB", "TMB", "rstan", "StanHeaders", "tmbstan", "Rcpp",
            "RcppEigen", "RcppParallel", "BH", "inline", "Matrix")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "-")
  cat(sprintf("%-12s: %s\n", p, v))
}
cat("stan version: ", rstan::stan_version(), "\n", sep = "")

set.seed(7)
y <- rnorm(20, mean = 3, sd = 1)

mk <- function(with_prior) {
  f <- function(pars) {
    mu <- pars$mu
    nll <- -sum(RTMB::dnorm(y, mu, 1, log = TRUE))
    if (with_prior) nll <- nll - RTMB::dnorm(mu, 0, 0.05, log = TRUE)
    nll
  }
  o <- RTMB::MakeADFun(f, list(mu = 0), silent = TRUE)
  o$fn(o$par)          # populate last.par.best
  o
}

o_flat <- mk(FALSE)
o_prior <- mk(TRUE)

cat("\nfn at mu = 0 : flat=", o_flat$fn(0), "  prior=", o_prior$fn(0), "\n",
    sep = "")
cat("fn at mu = 3 : flat=", o_flat$fn(3), "  prior=", o_prior$fn(3), "\n",
    sep = "")
cat("gr at mu = 1 : flat=", o_flat$gr(1), "  prior=", o_prior$gr(1), "\n",
    sep = "")
cat("DLL          : flat=", o_flat$env$DLL, " prior=", o_prior$env$DLL, "\n",
    sep = "")
cat("ADFun ptr    : flat=", format(o_flat$env$ADFun$ptr),
    " prior=", format(o_prior$env$ADFun$ptr), "\n", sep = "")
cat("ptr identical: ", identical(o_flat$env$ADFun$ptr, o_prior$env$ADFun$ptr),
    "\n", sep = "")

# analytic truth: 20 observations, sd 1, so the flat posterior is
# N(ybar, 1/sqrt(20)) = sd 0.2236; with the N(0, 0.05) prior the posterior sd
# is 1/sqrt(20 + 1/0.05^2) = 0.0471 and the mean is pulled to about 0.06.
cat("\nanalytic flat  : mean=", mean(y), " sd=", 1 / sqrt(20), "\n", sep = "")
post_v <- 1 / (20 + 1 / 0.05^2)
cat("analytic prior : mean=", post_v * sum(y), " sd=", sqrt(post_v), "\n",
    sep = "")

run <- function(o, dbg) {
  s <- tmbstan::tmbstan(o, chains = 1, iter = 1000, warmup = 500, seed = 42,
                        init = "0", debug = dbg,
                        refresh = 0, control = list(adapt_delta = 0.9))
  a <- rstan::extract(s, permuted = FALSE)
  cat("stan par names: ", paste(dimnames(a)[[3]], collapse = " "), "\n",
      sep = "")
  as.numeric(a[, 1, 1])
}

for (dbg in c(FALSE, TRUE)) {
  cat("\n### debug = ", dbg, " (", if (dbg) "R closures" else "tape pointer",
      ")\n", sep = "")
  a <- run(o_flat, dbg)
  b <- run(o_prior, dbg)
  cat("flat : mean=", mean(a), " sd=", sd(a), "\n", sep = "")
  cat("prior: mean=", mean(b), " sd=", sd(b), "\n", sep = "")
  cat("max|flat - prior| = ", max(abs(a - b)), "   identical = ",
      identical(a, b), "\n", sep = "")
}

cat("\n##### done\n")
