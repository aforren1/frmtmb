# Punch round 1: brms is the tiebreaker for where an ordinal fit's
# thresholds and cs() coefficients sit in fixef(), vcov() and
# summary()$fixed. MAJOR 1's table covers `ord ~ x`; the lane's own
# fixture also has `cs(z)`, and nothing measured said where `bcs_` goes.
#
# One R process, brms only, no frmtmb on the path.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"

set.seed(20260917)
n <- 120
dd <- data.frame(x = rnorm(n), z = rnorm(n),
                 g = factor(rep(1:12, each = 10)))
dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n),
                     c(-Inf, -0.3, 0.8, Inf), labels = 1:3),
                 ordered = TRUE)

ctl <- list(chains = 2, iter = 1000, refresh = 0, seed = 20260918,
            silent = 2, backend = "rstan", data = dd)
out <- list()
for (nm in c("cs", "plain")) {
  f <- if (nm == "cs") bf(ord ~ x + cs(z)) else bf(ord ~ x)
  fam <- if (nm == "cs") sratio() else cumulative()
  fit <- do.call(brm, c(list(formula = f, family = fam), ctl))
  s <- summary(fit)
  out[[nm]] <- list(
    fixef_rows = rownames(fixef(fit)),
    fixef_cols = colnames(fixef(fit)),
    vcov_dimnm = dimnames(vcov(fit)),
    fixed_rows = rownames(s$fixed),
    spec_pars = s$spec_pars,
    variables = variables(fit),
    fitted_dim = dim(fitted(fit)),
    fitted_dimnm = dimnames(fitted(fit)),
    fitted_range = range(fitted(fit)),
    fitted_q = range(fitted(fit)[, 3:4, ]),
    predict_dim = dim(predict(fit)),
    predict_dimnm = dimnames(predict(fit))
  )
  saveRDS(out, file.path(TREE, "dev/shapes-p1-brmsord.rds"))
  cat("== done", nm, "\n"); flush.console()
}
str(out)
