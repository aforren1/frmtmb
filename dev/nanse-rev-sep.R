# Reviewer: separated logistic regressions (test-predfix's design and a
# plain one). What the user sees: every warning, the SEs, glm's SEs.
#   Rscript dev/nanse-rev-sep.R base|lane|merge
args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
show <- function(lab, r) {
  se <- suppressWarnings(sqrt(diag(vcov(r$v))))
  cat(sprintf("%s: code %d logLik %.3g coef %s SE %s\n", lab,
              r$v$opt$convergence, as.numeric(logLik(r$v)),
              paste(signif(fixef(r$v, flatten = TRUE), 3), collapse = " "),
              paste(signif(se, 3), collapse = " ")))
  for (x in r$w) cat("    warn:", substr(x, 1, 140), "\n")
}
for (seed in 511:514) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
  d$yb <- as.integer(d$z > 0)
  show(paste("predfix seed", seed, "default"),
       cap(frm(yb ~ z + x, family = bernoulli(), data = d)))
}
for (seed in 1:5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(80), z = rnorm(80))
  d$y <- as.integer(d$x + 0.3 * d$z > 0)
  g <- suppressWarnings(glm(y ~ x + z, binomial, data = d))
  cat(sprintf("glm seed %d: coef %s SE %s\n", seed,
              paste(signif(coef(g), 3), collapse = " "),
              paste(signif(sqrt(diag(vcov(g))), 3), collapse = " ")))
  show(paste("plain separated seed", seed),
       cap(frm(y ~ x + z, family = bernoulli(), data = d)))
}
