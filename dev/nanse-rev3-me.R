# Reviewer, punch round 2: interop_vcov() hands marginaleffects the
# displayed covariance (lost rows NaN). What does marginaleffects report
# along a lost direction, and for a quantity that does not move along
# one? Fit 1: predfix seed 514 (separated, code 0, every coefficient
# lost). Fit 2: y ~ x + z with z held by a bound (ub), a lost
# coefficient that is not a flat direction. Fit 3: a mixed model whose
# (1 + x | g) block is flat (test-portability's fixture), where the
# fixed effects are determined.
#   Rscript dev/nanse-rev3-me.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(marginaleffects)})
cat("arm", arm, "from", find.package("frmtmb"), "marginaleffects",
    as.character(packageVersion("marginaleffects")), "\n")
ns <- asNamespace("frmtmb")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "err")
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
show <- function(lab, r) {
  if (inherits(r$v, "err")) {
    cat(sprintf("  %-34s ERROR %s\n", lab, substr(r$v, 1, 120)))
    return(invisible())
  }
  se <- r$v$std.error
  cat(sprintf("  %-34s std.error %s | finite %d of %d | warnings %d%s\n",
              lab, paste(signif(head(se, 3), 3), collapse = " "),
              sum(is.finite(se)), length(se), length(r$w),
              if (length(r$w)) paste0(" [", substr(r$w[1], 1, 60), "]") else ""))
}
set.seed(514)
d1 <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
d1$yb <- as.integer(d1$z > 0)
f1 <- cap(frm(yb ~ z + x, family = bernoulli(), data = d1))$v
cat("fit 1 lost:", paste(names(ns$sdr_of(f1)$se_lost), collapse = ","), "\n")
show("predictions(type = link)", cap(predictions(f1, type = "link")))
show("avg_slopes(variables = z)", cap(avg_slopes(f1, variables = "z")))
set.seed(101)
d2 <- data.frame(x = rnorm(300), z = rnorm(300))
d2$y <- rnorm(300, 1 + 0.5 * d2$x + 2 * d2$z, 1)
f2 <- suppressWarnings(frm(y ~ x + z, data = d2,
                           prior = set_prior("", class = "b", coef = "z",
                                             ub = 0.1)))
cat("fit 2 lost:", paste(sprintf("%s(%s)", names(ns$sdr_of(f2)$se_lost),
                                 ns$sdr_of(f2)$se_lost), collapse = ","), "\n")
show("predictions()", cap(predictions(f2)))
show("avg_slopes(variables = x)", cap(avg_slopes(f2, variables = "x")))
show("frm_linpred(se.fit)", list(v = list(std.error = cap(frm_linpred(
  f2, newdata = d2[1:3, ], se.fit = TRUE))$v$se.fit), w = character()))
set.seed(1)
n <- 60
d3 <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(n / 10, 10), y = 0)
d3$y <- frm_simulate(bf(y ~ x) + gaussian(), d3,
                     newparams = list(b_Intercept = 1, b_x = 2, sigma = 1),
                     nsim = 1, seed = 1001L)[[1]]
f3 <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)), data = d3))
cat("fit 3 lost:", paste(names(ns$sdr_of(f3)$se_lost), collapse = ","), "\n")
show("predictions(re.form = NA)", cap(predictions(f3, re.form = NA)))
show("avg_slopes(variables = x)", cap(avg_slopes(f3, variables = "x")))
