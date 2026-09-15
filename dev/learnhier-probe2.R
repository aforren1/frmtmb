# Lane `learnhier`: the rlddm block's row order, the confint_varcorr()
# row names, and what a pinned `bias` leaves in the parameter vector.
# Toy sizes, so it costs seconds.
source("dev/learnhier-env.R")
source("dev/learnhier-sim.R")
suppressPackageStartupMessages(library(frmtmb.learn))

cat("eigenvalues of the rlddm truth correlation:\n")
print(round(eigen(lh_rlddm_truth$R, only.values = TRUE)$values, 4))
cat("eigenvalues of the bandit truth correlation:\n")
print(round(eigen(lh_bandit_truth$R, only.values = TRUE)$values, 4))

d <- lh_rlddm_data(4001L, ns = 10L, nt = 80L)
cat("\nrows:", nrow(d), " floors mean:", mean(attr(d, "own_floor")), "\n")
cat("induced block stats (no fit):\n")
print(round(lh_block_stats(attr(d, "dev_fitted")), 4))
cat("drawn block stats (no fit):\n")
print(round(lh_block_stats(attr(d, "dev_drawn")), 4))

form <- frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) +
                     ndt_group(id) ~ 1 + (1 | p | id),
                   drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
                   ndt ~ 1 + (1 | p | id), bias = 0.5)
fit <- frmtmb::frm(form, family = rlddm(subject = id, trial = trial),
                   data = d)
cat("\nVarCorr rownames:\n")
print(rownames(frmtmb::VarCorr(fit)[[1L]]))
cat("\nfixef:\n")
print(unlist(frmtmb::fixef(fit)))
cat("\npar name table:\n")
print(table(names(fit$obj$env$last.par.best)))
cat("\nconfint rownames:\n")
print(rownames(suppressWarnings(stats::confint(fit))))
cat("\nconfint_varcorr:\n")
print(suppressWarnings(frmtmb::confint_varcorr(fit)))
cat("\nndt_bound floors:\n")
bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
print(utils::str(bd, max.level = 1L))
cat("\nb entries, first 12 (level-major, 4 per level):\n")
p <- fit$obj$env$last.par.best
print(round(unname(p[names(p) == "b"])[1:12], 4))
cat("\nranef head:\n")
print(utils::head(frmtmb::ranef(fit)[[1L]], 3L))
