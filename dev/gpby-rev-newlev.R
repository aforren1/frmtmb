# Reviewer: every frmtmb.spline route whose allow_new_levels refusal the
# lane removed, at an unseen level of (1 + t | g), against closed forms
# built from VarCorr() and vcov() and a Monte Carlo of the curve, its
# slope and its zero crossing. Lane build only (base refuses).
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
set.seed(23)
ng <- 25
d <- data.frame(g = factor(rep(seq_len(ng), each = 12)),
                t = rep(seq(0, 1, length.out = 12), ng))
u0 <- rnorm(ng, 0, 0.5); u1 <- rnorm(ng, 0, 0.3)
d$y <- -1 + 2 * d$t + u0[d$g] + u1[d$g] * d$t + rnorm(nrow(d), 0, 0.2)
fit <- frm(bf(y ~ t + (1 + t | g)), data = d)
vc <- varcorr_matrices <- NULL
Sg <- frmtmb::VarCorr(fit)$g$cov
Sg <- Sg[, "Estimate", ]
cat("Sigma_g:\n"); print(Sg)
Vb <- vcov(fit)[c("Intercept", "t"), c("Intercept", "t")]
if (anyNA(Vb)) Vb <- vcov(fit)[1:2, 1:2]
be <- fixef(fit)[1:2, "Estimate"]
tt <- seq(0, 1, length.out = 21)
lev <- c(levels(d$g), "n1", "n2")
nd1 <- data.frame(t = tt, g = factor("n1", levels = lev))
nd2 <- data.frame(t = tt, g = factor("n2", levels = lev))
X <- cbind(1, tt)
# (a) the band
cv <- frm_curve(fit, newdata = nd1, re_formula = NULL,
                allow_new_levels = TRUE, nsim = 20000, seed = 1)
Scf <- X %*% Vb %*% t(X) + X %*% Sg %*% t(X)
cat(sprintf("BAND Sigma vs closed form max rel %.3e | crit %.4f\n",
            max(abs(attr(cv, "Sigma") - Scf)) / max(abs(Scf)),
            cv$.crit_sim[1]))
set.seed(99)
nmc <- 20000
bd <- MASS::mvrnorm(nmc, be, Vb)
ud <- MASS::mvrnorm(nmc, c(0, 0), Sg)
curves <- (bd + ud) %*% t(X)
inband <- apply(curves, 1, function(r) all(r >= cv$.lower_sim &
                                            r <= cv$.upper_sim))
cat(sprintf("BAND MC whole-curve coverage of draws %.4f (nominal 0.95, mcse %.4f)\n",
            mean(inband), sqrt(0.95 * 0.05 / nmc)))
# (b) the derivative
dv <- frm_curve_deriv(fit, var = "t", newdata = nd1, re_formula = NULL,
                      allow_new_levels = TRUE, simultaneous = FALSE)
cat(sprintf("DERIV se %.6f vs closed form sqrt(Var(b_t) + Sg[2,2]) %.6f, ratio %.8f | seen-level-free se %.6f\n",
            dv$.se[5], sqrt(Vb[2, 2] + Sg[2, 2]),
            dv$.se[5] / sqrt(Vb[2, 2] + Sg[2, 2]), sqrt(Vb[2, 2])))
# (c) the crossing of 0 at the unseen level
ft <- frm_curve_feature(fit, var = "t", type = "crossing", at = 0,
                        newdata = nd1, re_formula = NULL,
                        allow_new_levels = TRUE)
ts <- ft$.estimate[1]
x0 <- c(1, ts)
se_cf <- sqrt(drop(t(x0) %*% (Vb + Sg) %*% x0)) / abs(be[2])
roots <- -(bd[, 1] + ud[, 1]) / (bd[, 2] + ud[, 2])
ok <- roots > -1 & roots < 2
cat(sprintf(paste0("CROSS root %.5f (analytic %.5f) | se %.6f | delta ",
                   "closed form %.6f ratio %.8f | MC sd %.6f, MAD %.6f ",
                   "(%.3f of draws in (-1, 2)) | value_se %.6f vs %.6f\n"),
            ts, -be[1] / be[2], ft$.se[1], se_cf, ft$.se[1] / se_cf,
            sd(roots[ok]), mad(roots[ok]), mean(ok), ft$.value_se[1],
            sqrt(drop(t(x0) %*% (Vb + Sg) %*% x0))))
# (d) differences
c12 <- frm_curve(fit, newdata = nd1, contrast = nd2, re_formula = NULL,
                 allow_new_levels = TRUE, simultaneous = FALSE)
c11 <- frm_curve(fit, newdata = nd1, contrast = nd1, re_formula = NULL,
                 allow_new_levels = TRUE, simultaneous = FALSE)
cat(sprintf(paste0("DIFF n1 - n2 se vs sqrt(2 x' Sg x) max rel %.3e | ",
                   "n1 - n1 max se %.3e\n"),
            max(abs(c12$.se / sqrt(2 * rowSums((X %*% Sg) * X)) - 1)),
            max(c11$.se)))
dd <- frm_curve_deriv(fit, var = "t", newdata = nd1, contrast = nd2,
                      re_formula = NULL, allow_new_levels = TRUE,
                      simultaneous = FALSE)
cat(sprintf("DIFF deriv n1 - n2 se %.6f vs sqrt(2 Sg[2,2]) %.6f\n",
            dd$.se[5], sqrt(2 * Sg[2, 2])))
# (e) an exact gp() past the data: zero crossing and slope against a
# Monte Carlo of the field drawn from the fitted model's joint law
set.seed(5)
x <- sort(runif(50, 0, 4))
dg <- data.frame(x = x, y = sin(1.3 * x) + rnorm(50, 0, 0.15))
fg <- frm(bf(y ~ gp(x)), data = dg)
grid <- data.frame(x = seq(3.5, 5.5, length.out = 41))
lb <- frm_lp_basis(fg, newdata = grid, extra_cov = TRUE)
A <- as.matrix(lb$A)
S <- A %*% lb$V %*% t(A) + lb$extra_cov
cvg <- frm_curve(fg, newdata = grid, nsim = 20000, seed = 2)
cat(sprintf("GP band Sigma == A V A' + extra_cov max rel %.3e\n",
            max(abs(attr(cvg, "Sigma") - S)) / max(abs(S))))
ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))
set.seed(7)
dr <- lb$eta + L %*% matrix(rnorm(41 * nmc), 41, nmc)
inb <- colMeans(dr >= cvg$.lower_sim & dr <= cvg$.upper_sim) == 1
cat(sprintf("GP band MC whole-curve coverage of draws %.4f\n", mean(inb)))
fgx <- frm_curve_feature(fg, var = "x", type = "crossing", at = 0,
                         newdata = grid)
cat("GP crossings:", paste(sprintf("%.4f se %.4f", fgx$.estimate, fgx$.se),
                           collapse = "; "), "\n")
DONE <- TRUE
cat("DONE\n")
