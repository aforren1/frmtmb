# Reviewer, punch round 2: how often does the noise term of se_tier3()
# take the SEs of identified fixed effects? Design of
# dev/nanse-rev3-falseloss.R (gaussian, y ~ x + f + (1 + x | g2), 40-level
# f, no g2 variation), seeds 1 to 20, plus a smaller one without the
# factor (y ~ x + (1 + x | g2), 120 rows). Per seed: fixed effects lost on
# the lane, and lme4's SE for the intercept for comparison.
#   Rscript dev/rel068-falseloss-sweep.R release (the reviewer's script with a
#   release arm for the 0.68.0 merge; lane|merge as before)
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  release = c("C:/Users/adf44/source/r/rellib-r6"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
tab <- NULL
for (design in c("f40", "plain")) {
  for (s in 1:20) {
    set.seed(s)
    n <- if (design == "f40") 480 else 120
    d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                    g2 = factor(rep(1:6, length.out = n)))
    d$y <- 1 + 0.5 * d$x + (if (design == "f40") rnorm(40, 0, 0.5)[d$f]
                            else 0) + rnorm(n)
    fo <- if (design == "f40") y ~ x + f + (1 + x | g2) else
      y ~ x + (1 + x | g2)
    w <- character()
    fit <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
    lost <- suppressWarnings(ns$sdr_of(fit))$se_lost
    om <- ns$outer_par_map(fit)
    beta_nm <- om$names[om$comp == "beta"]
    nb <- sum(names(lost) %in% beta_nm)
    se_i <- suppressWarnings(sqrt(diag(vcov(fit))))[1]
    m <- suppressMessages(suppressWarnings(
      lme4::lmer(fo, data = d, REML = FALSE)))
    tab <- rbind(tab, data.frame(design, seed = s, code = fit$opt$convergence,
                                 n_lost = length(lost), beta_lost = nb,
                                 se_int = se_i,
                                 lme4_se_int = sqrt(vcov(m)[1, 1]),
                                 singular = lme4::isSingular(m),
                                 warned = any(grepl("Standard errors are not",
                                                    w))))
  }
}
print(tab, row.names = FALSE)
cat("\nfits with identified fixed effects lost:\n")
print(aggregate(cbind(fits = 1, beta_lost_fits = beta_lost > 0,
                      singular = singular) ~ design, tab, sum))
