# Punch round 2: brms is the tiebreaker for three shapes.
#   (D) coef() on a mixed ordinal fit, for a thres_minus_eta family
#       (cumulative) and an eta_minus_thres one (acat);
#   (C) the joint draws of a rescor multivariate gaussian;
#   (minor) posterior_predict/predict(summary = FALSE) dimnames and
#       ranef() column names.
# One R process, brms only.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
set.seed(20260921)
n <- 150
dd <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
u <- rnorm(15, 0, 0.8)
dd$ord <- factor(cut(0.8 * dd$x + u[dd$g] + rlogis(n), c(-Inf, -0.5, 0.8, Inf),
                     labels = 1:3), ordered = TRUE)
e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.7, 0.7, 1), 2))
dd$y1 <- 1 + 0.5 * dd$x + e[, 1]
dd$y2 <- -0.3 * dd$x + e[, 2]
rownames(dd) <- paste0("case", seq_len(n))
ctl <- list(chains = 2, iter = 1000, refresh = 0, seed = 20260921,
            silent = 2, backend = "rstan", data = dd)
out <- list()
for (fam in c("cumulative", "acat")) {
  fit <- do.call(brm, c(list(formula = bf(ord ~ x + (1 | g)),
                             family = get(fam)()), ctl))
  cf <- coef(fit, summary = FALSE)$g
  fe <- fixef(fit, summary = FALSE)
  re <- ranef(fit, summary = FALSE)$g
  out[[fam]] <- list(
    coef_names = dimnames(coef(fit)$g)[[3]],
    ranef_names = dimnames(ranef(fit)$g)[[3]],
    # which formula reproduces brms's own coef, draw by draw
    tau_minus_r = max(abs(cf[, 1, "Intercept[1]"] -
                            (fe[, "Intercept[1]"] - re[, 1, "Intercept"]))),
    r_minus_tau = max(abs(cf[, 1, "Intercept[1]"] -
                            (re[, 1, "Intercept"] - fe[, "Intercept[1]"]))))
  cat("==", fam, "\n"); flush.console()
}
fit <- do.call(brm, c(list(formula = bf(mvbind(y1, y2) ~ x) +
                             set_rescor(TRUE), family = gaussian()), ctl))
pp <- posterior_predict(fit)
pr <- predict(fit, summary = FALSE)
out$mv <- list(
  rescor = summary(fit)$rescor_pars,
  pp_dim = dim(pp), pp_dimnames = dimnames(pp),
  pr_dim = dim(pr), pr_dimnames = dimnames(pr),
  draw_cor = mean(vapply(seq_len(n), function(i) cor(pp[, i, 1], pp[, i, 2]),
                         0)))
fit1 <- do.call(brm, c(list(formula = bf(y1 ~ x + (1 | g))), ctl))
out$uni <- list(pp_dimnames = dimnames(posterior_predict(fit1)),
                pr_dimnames = dimnames(predict(fit1, summary = FALSE)),
                ranef_names = dimnames(ranef(fit1)$g)[[3]],
                coef_names = dimnames(coef(fit1)$g)[[3]])
saveRDS(out, file.path(TREE, "dev/shapes-p2-brmsref.rds"))
str(out)
