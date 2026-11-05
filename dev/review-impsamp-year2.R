# B1's other spellings: an uncentered column in a dispersion's regression, and in
# a random slope. The centered column is the reference: the model is the same
# (us() is invariant under an affine map of the column), so logLik must agree.
# Usage: Rscript dev/review-impsamp-year2.R after|mutant
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- commandArgs(TRUE)[[1]]
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
cat("lib", find.package("frmtmb"), "\n")
cases <- list(
  shape = function(s, col) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10))
    z <- rnorm(150)
    y <- rnbinom(150, size = exp(0.5 + 0.4 * z), mu = exp(1 + rnorm(15, 0, 0.6)[g]))
    d <- data.frame(y, g, yr = round(2000 + 5 * z))
    d$w <- if (col == "raw") d$yr else d$yr - mean(d$yr)
    list(f = bf(y ~ 1 + (1 | g), shape ~ w), fam = negbinomial(), d = d)
  },
  slope = function(s, col) {
    set.seed(s)
    g <- factor(rep(seq_len(25L), each = 8))
    z <- rnorm(200)
    b <- rnorm(25, 0, 0.3)
    y <- rpois(200, exp(0.5 + 0.3 * z + rnorm(25, 0, 0.6)[g] + b[g] * z))
    d <- data.frame(y, g, yr = 2000 + 5 * z)
    d$w <- if (col == "raw") d$yr else d$yr - mean(d$yr)
    list(f = bf(y ~ w + (1 + w | g)), fam = poisson(), d = d)
  })
for (cn in names(cases)) {
  for (s in 1:4) {
    r <- list()
    for (col in c("raw", "centered")) {
      de <- cases[[cn]](s, col)
      w <- character(0)
      t0 <- Sys.time()
      fit <- withCallingHandlers(
        tryCatch(frm(de$f, family = de$fam, data = de$d, importance = 300L),
                 error = function(e) e),
        warning = function(x) {
          w <<- c(w, substr(conditionMessage(x), 1, 60))
          invokeRestart("muffleWarning")
        })
      lap <- suppressWarnings(tryCatch(frm(de$f, family = de$fam, data = de$d),
                                       error = function(e) NULL))
      r[[col]] <- if (inherits(fit, "error")) {
        sprintf("STOP %s", substr(conditionMessage(fit), 1, 70))
      } else {
        se <- tryCatch(sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed)), error = function(e) NA)
        sprintf("ll %.3f (lap %.3f) rounds %d grad %.3g ess %.2f se_slope %.4g warns %d %s",
                as.numeric(logLik(fit)), if (is.null(lap)) NA else as.numeric(logLik(lap)),
                fit$importance$rounds, fit$importance$grad, fit$importance$ess_min,
                se[[2]], length(w), paste(w, collapse = ";"))
      }
    }
    cat(cn, s, "raw     :", r$raw, "\n")
    cat(cn, s, "centered:", r$centered, "\n")
  }
}
