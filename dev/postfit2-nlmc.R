# Lane postfit2: the nonlinear Wald band of conditional_effects()
# against Monte Carlo over the fit's covariance.
#
#   Rscript dev/postfit2-nlmc.R > dev/postfit2-log/nlmc.txt
#
# Three checks, each on the model of dev/postfit2-brms-compare.R B4
# (data seed 1, n = 150, y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1):
# 1. The Jacobian frm_lp_basis() tapes against central finite
#    differences of frm_linpred() in each coefficient.
# 2. se__ against the SD of R curves at coefficients drawn from
#    N(chat, V), two seeds, R = 20000 each. A curved body makes the
#    true SD differ from the first-order one, so a gap that survives a
#    larger R is curvature, not Monte Carlo error.
# 3. The same with V / 100: curvature shrinks with the spread and the
#    delta method becomes exact, so here se__ / 10 must match the MC SD
#    to within its Monte Carlo error. This is the check of the
#    implementation; 2 is the size of the approximation.
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(1)
n <- 150
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(2 * pi * d$x) + d$z * (d$f == "b") + rnorm(n, 0, 0.3)
d$y2 <- 2 * exp(0.8 * d$z) * (1 + 0.3 * (d$f == "b")) + rnorm(n, 0, 0.3)
fit <- frm(bf(y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1, nl = TRUE),
           family = gaussian(), data = d)
ce <- conditional_effects(fit, "z:f", resolution = 50)[[1]]
nd <- ce
lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
pm <- frmtmb:::joint_pos_map(frmtmb:::get_joint_cov(fit))
chat <- vapply(lb$coef_pos, function(k) {
  fit$estimates[[pm$comp[k]]][pm$idx[k]]
}, 0)
eta_at <- function(cv) {
  est <- fit$estimates
  for (j in seq_along(lb$coef_pos)) {
    k <- lb$coef_pos[j]
    est[[pm$comp[k]]][pm$idx[k]] <- cv[j]
  }
  fz <- fit
  fz$estimates <- est
  as.vector(frm_linpred(fz, newdata = nd, type = "link", dpar = "mu",
                        re_formula = NA))
}
say("coefficients: %s", paste(lb$coef_names, collapse = ", "))
say("se__ equals sqrt(diag(A V A')): max rel diff %.3g",
    max(abs(ce$se__ - sqrt(rowSums((lb$A %*% lb$V) * lb$A))) / ce$se__))

## 1. Jacobian against finite differences
h <- 1e-5
fd <- sapply(seq_along(chat), function(j) {
  e <- replace(numeric(length(chat)), j, h)
  (eta_at(chat + e) - eta_at(chat - e)) / (2 * h)
})
say("1. Jacobian vs central differences (h = %g): max |diff| / max |A| = %.3g",
    h, max(abs(fd - lb$A)) / max(abs(lb$A)))

## 2 and 3. Monte Carlo
mc <- function(V, R, seed) {
  set.seed(seed)
  L <- t(chol(V))
  sims <- matrix(NA_real_, R, nrow(nd))
  for (r in seq_len(R)) {
    sims[r, ] <- eta_at(chat + as.vector(L %*% rnorm(length(chat))))
  }
  apply(sims, 2, stats::sd)
}
R <- 20000
for (scale in c(1, 0.01)) {
  for (seed in c(20260929, 20260930)) {
    s <- mc(lb$V * scale, R, seed)
    se <- ce$se__ * sqrt(scale)
    ratio <- se / s
    z <- (se - s) / (s / sqrt(2 * (R - 1)))
    say("%s V = %s, R = %d, seed %d: se/MCsd min %.4f median %.4f max %.4f; |z| max %.2f, mean z %.2f; MC error of one SD %.4f",
        if (scale == 1) "2." else "3.", format(scale), R, seed,
        min(ratio), median(ratio), max(ratio), max(abs(z)), mean(z),
        1 / sqrt(2 * (R - 1)))
  }
}
