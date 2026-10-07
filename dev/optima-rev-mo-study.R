# Reviewer of lane optima, claim 1: the headline count, recomputed
# independently of dev/optima-mo-study.R. brms_monotonic's data code,
# `ls ~ mo(income) * age`, seeds 1 to 200, against the exact maximum
# by enumeration over signs and supports (test-optima.R's construction).
#   Rscript dev/optima-rev-mo-study.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
dat <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
exact <- function(d) {
  code <- as.integer(d$income) - 1L
  S <- sapply(1:3, function(k) as.numeric(code >= k))
  n <- nrow(d)
  best <- -Inf
  for (s1 in c(-1, 1)) for (s2 in c(-1, 1)) for (m1 in 0:7) for (m2 in 0:7) {
    a1 <- which(bitwAnd(m1, c(1L, 2L, 4L)) > 0)
    a2 <- which(bitwAnd(m2, c(1L, 2L, 4L)) > 0)
    X <- cbind(1, d$age, S[, a1, drop = FALSE],
               (S * d$age)[, a2, drop = FALSE])
    f <- stats::lm.fit(X, d$ls)
    cf <- f$coefficients[-(1:2)]
    if (anyNA(cf)) next
    sg <- c(rep(s1, length(a1)), rep(s2, length(a2)))
    if (any(sg * cf < 0)) next
    best <- max(best, -(n / 2) * (log(2 * pi * sum(f$residuals^2) / n) + 1))
  }
  best
}
gap <- code <- ev <- nw <- nse <- numeric(200)
for (s in 1:200) {
  d <- dat(s)
  w <- character()
  f <- withCallingHandlers(frm(ls ~ mo(income) * age, data = d),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  gap[s] <- exact(d) - as.numeric(logLik(f))
  code[s] <- f$opt$convergence
  ev[s] <- f$opt$evals
  nw[s] <- length(w)
  nse[s] <- sum(!is.finite(fixef(f)[, "Est.Error"]))
}
cat("gap <= 1e-6:", sum(gap <= 1e-6), "; > 1e-2:", sum(gap > 1e-2),
    "; max", format(max(gap), digits = 4), "; min", format(min(gap), digits = 3),
    "\ncodes != 0:", sum(code != 0), "; fits with a warning:", sum(nw > 0),
    "; fits with a non-finite fixef SE:", sum(nse > 0),
    "; evaluations", sum(ev), "\n")
cat("seeds above 1e-6:", which(gap > 1e-6)[seq_len(min(30, sum(gap > 1e-6)))],
    "\n")
