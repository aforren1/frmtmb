# Reviewer, punch round 2: on dev/nanse-rev3-deferred-paths.R's natural
# deferral fixture (seed 11), the lane loses (Intercept) and f2..f40 as
# well as the flat g2 block. Is that a false loss? Reference: lme4 (ML)
# on the same model, and base's sdreport(). Also the tier-3 numbers:
# the noise norm of the finite-difference Hessian, flat_thr and tau.
#   Rscript dev/nanse-rev3-falseloss.R lane|base [seed]
args <- commandArgs(TRUE)
arm <- args[1]
s <- if (length(args) > 1) as.integer(args[2]) else 11L
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "seed", s, "\n")
ns <- asNamespace("frmtmb")
set.seed(s)
n <- 480
d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                g2 = factor(rep(1:6, length.out = n)))
d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$f] + rnorm(n)
w <- character()
fit <- withCallingHandlers(frm(y ~ x + f + (1 + x | g2), data = d),
  warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
V <- suppressWarnings(vcov(fit, full = TRUE))
se <- suppressWarnings(sqrt(diag(V)))
nm <- ns$outer_par_names(fit)
names(se) <- nm
lost <- ns$sdr_of(fit)$se_lost
cat("code", fit$opt$convergence, "| lost", length(lost), ":",
    paste(head(names(lost), 6), collapse = " "), "...\n")
cat("theta:", signif(fit$opt$par[grepl("^theta", nm)], 4), "\n")
cat("SE (Intercept) x f2 f40 sigma:",
    signif(se[c("(Intercept)", "x", "f2", "f40", "sigma_(Intercept)")], 5),
    "\n")
h <- ns$fit_outer_hessian(fit)
H <- h$H
E <- h$E
D <- sqrt(abs(diag(H)))
S <- H / outer(D, D)
ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
cat("unit-diag eigenvalues: smallest 6", signif(tail(ev, 6), 3),
    "largest", signif(ev[1], 3), "\n")
Es <- E / outer(D, D)
cat("noise ||E||_F (unit-diag) =", signif(sqrt(sum(Es^2)), 3),
    "-> flat_thr = 10 * that =", signif(10 * sqrt(sum(Es^2)), 3), "\n")
if (requireNamespace("lme4", quietly = TRUE)) {
  m <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ x + f + (1 + x | g2), data = d, REML = FALSE)))
  sm <- sqrt(diag(as.matrix(vcov(m))))
  cat("lme4 ML: logLik", format(logLik(m), digits = 10), "frmtmb",
      format(logLik(fit), digits = 10), "\n")
  cat("lme4 SE (Intercept) x f2 f40:",
      signif(sm[c("(Intercept)", "x", "f2", "f40")], 5), "singular",
      lme4::isSingular(m), "\n")
}
# where the noise sits: per parameter, its unit-diagonal noise row norm
rn <- sqrt(rowSums(Es^2))
o <- order(-rn)[1:5]
cat("largest noise rows:", paste(sprintf("%s=%.3g (|H_ii| %.3g, row max %.3g)",
    nm[o], rn[o], abs(diag(H))[o], apply(abs(H), 1, max)[o]), collapse = "; "),
    "\n")
cat("noise norm without those rows:",
    signif(sqrt(sum(Es[-o[1], -o[1]]^2)), 3), "\n")
