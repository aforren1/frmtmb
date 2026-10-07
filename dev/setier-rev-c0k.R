# Reviewer of lane setier: c0 + exp(a)^k with a ~ 1 + (1 | g) (lane nanse
# review's c0k model, seed 77). The Laplace likelihood is invariant to
# (a0, k, b) -> (s a0, k / s, s b), so a_(Intercept), k_(Intercept) and
# theta_1 (log sd of b) lie on a curved ridge. What does tier 3 remove,
# what does se_curvature_real() rescue, and is the objective constant
# along the ridge?
#   Rscript dev/setier-rev-c0k.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  merge = c("C:/Users/adf44/source/r/setier-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(77)
G <- 8
dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(G, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
f <- suppressMessages(suppressWarnings(
  frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE),
      data = dn)))
nm <- ns$outer_par_names(f)
p <- f$opt$par
names(p) <- nm
cat("par:", paste0(nm, "=", signif(p, 5), collapse = " "), "\n")
se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
cat("SE :", paste0(nm, "=", signif(se, 4), collapse = " "), "\n")
cat("lost:", paste(names(ns$sdr_of(f)$se_lost), ns$sdr_of(f)$se_lost),
    "\n")
f0 <- f$obj$fn(p)
ia <- match("a_(Intercept)", nm)
ik <- match("k_(Intercept)", nm)
it <- match("theta_1", nm)
for (s in c(0.5, 0.9, 1.1, 2)) {
  q <- p
  q[ia] <- p[ia] * s
  q[ik] <- p[ik] / s
  q[it] <- p[it] + log(s)
  cat(sprintf("ridge s = %.1f: nll change %.3g\n", s, f$obj$fn(q) - f0))
}
h <- ns$fit_outer_hessian(f)
H <- h$H
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("unit-diag ev:", signif(e$values, 3), "\n")
cat("max:", signif(max(abs(e$values)), 3), " 1e-9 * max:",
    signif(1e-9 * max(abs(e$values)), 3), "\n")
Es <- h$E / outer(D, D)
nk <- sqrt(colSums((Es %*% e$vectors)^2))
cat("own noise ||Es v||:", signif(nk, 3), "\n")
for (k in which(e$values <= pmax(1e-9 * max(abs(e$values)), 3 * nk))) {
  cat("direction", k, "ev", signif(e$values[k], 3), "loads:",
      paste0(nm, "=", round(e$vectors[, k], 3), collapse = " "), "\n")
  if (exists("se_curvature_real", ns) && e$values[k] > 0) {
    dir <- e$vectors[, k] / D
    lam <- e$values[k]
    tol <- f$control$grad_tol %||% 1e-3
    d <- dir * sqrt(4 * tol / lam)
    cat("   probe step (natural units):", signif(d, 3), "\n")
    cat("   nll change +step", signif(f$obj$fn(p + d) - f0, 4), " -step",
        signif(f$obj$fn(p - d) - f0, 4), "\n")
    cat("   se_curvature_real:", ns$se_curvature_real(f, p, seq_along(p),
                                                       dir, lam), "\n")
  }
  # the exact ridge tangent at p, unit-diag coordinates
  tg <- numeric(length(p))
  tg[ia] <- p[ia]; tg[ik] <- -p[ik]; tg[it] <- 1
  tg <- tg * D
  tg <- tg / sqrt(sum(tg^2))
  cat("   |cos| with the exact ridge tangent:",
      signif(abs(sum(tg * e$vectors[, k])), 4), "\n")
}
