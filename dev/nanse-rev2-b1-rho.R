# Reviewer, punch round 1: attacks on the projection rule of se_tier3().
# Each design has a ridge (so tier 3 runs) beside parameters the data
# DO determine; those must keep the SE of the same model written
# without the ridge (the reference).
#   (generated from nanse-rev2-b1.R with cor 1 - 1e-7 and 1 - 1e-9)
#   A: c + d ridge beside near-collinear a ~ 0 + x1, b ~ 0 + x2
#      (cor 0.999 and 0.99999); exact AD Hessian (no random effects).
#   B: the same with (1 | g) in c, so the finite-difference Hessian and
#      its noise term decide.
#   C: a two-dimensional flat subspace with mixed loadings,
#      y ~ a + b + c + e, a ~ 0 + f (30 levels), b ~ 0 + h (3 levels),
#      c ~ 1, e ~ 0 + x, x correlated with f; every a, b, c lost, e kept.
#   Rscript dev/nanse-rev2-b1.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
se_of <- function(fit) {
  V <- suppressWarnings(vcov(fit, full = TRUE))
  s <- suppressWarnings(sqrt(diag(V)))
  names(s) <- ns$outer_par_names(fit)
  s
}
show <- function(lab, r, ref, keep) {
  f <- r$v
  s <- se_of(f)
  sr <- se_of(ref)
  lost <- ns$sdr_of(f)$se_lost
  rel <- s[keep] / sr[keep] - 1
  cat(sprintf("%s: code %d | lost %s | %s | max |SE/ref - 1| %s\n", lab,
              f$opt$convergence,
              if (length(lost)) paste(names(lost), collapse = ",") else "none",
              paste(sprintf("%s=%.4g (ref %.4g)", keep, s[keep], sr[keep]),
                    collapse = " "),
              format(max(abs(rel)), digits = 3)))
  for (x in r$w) cat("    warn:", substr(x, 1, 130), "\n")
}
for (rho in c(0.9999999, 0.999999999)) {
  set.seed(7)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rho * x1 + sqrt(1 - rho^2) * rnorm(n)
  g <- factor(rep(1:15, 20))
  d <- data.frame(x1, x2, g, y = 1 + 0.5 * x1 + 0.5 * x2 +
                    rnorm(15, 0, 0.5)[g] + rnorm(n))
  cat(sprintf("\n== cor(x1, x2) = %.6f\n", cor(x1, x2)))
  fA <- cap(frm(bf(y ~ a + b + c + dd, a ~ 0 + x1, b ~ 0 + x2, c ~ 1,
                   dd ~ 1, nl = TRUE), data = d))
  rA <- cap(frm(bf(y ~ a + b + c, a ~ 0 + x1, b ~ 0 + x2, c ~ 1, nl = TRUE),
                data = d))
  show("A exact", fA, rA$v, c("a_x1", "b_x2", "sigma_(Intercept)"))
  fB <- cap(frm(bf(y ~ a + b + c + dd, a ~ 0 + x1, b ~ 0 + x2,
                   c ~ 1 + (1 | g), dd ~ 1, nl = TRUE), data = d))
  rB <- cap(frm(bf(y ~ a + b + c, a ~ 0 + x1, b ~ 0 + x2, c ~ 1 + (1 | g),
                   nl = TRUE), data = d))
  show("B finite-difference", fB, rB$v,
       c("a_x1", "b_x2", "sigma_(Intercept)", "theta_1"))
}
