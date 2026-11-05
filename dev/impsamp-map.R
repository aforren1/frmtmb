# A constant dpar reaches the corrected objective through a map. Does the
# numeric path that sets the offsets and reads the ESS evaluate the same
# model as the tape? Compare, at the fit's own estimate, the sum over
# groups of the log-mean-exp from amat() with the tape's objective.
# Data seed 31. Usage: Rscript dev/impsamp-map.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(31)
g <- factor(rep(1:15, each = 8))
x <- rnorm(120)
y <- 0.3 + 0.5 * x + rnorm(15, 0, 0.7)[g] + 0.5 * rt(120, 5)
d <- data.frame(y = y, x = x, g = g)
fit <- suppressWarnings(frm(bf(y ~ x + (1 | g), nu = 5), family = student(),
                            data = d, importance = 200))
cat("outer par:", paste(names(fit$opt$par), collapse = " "), "\n")
cat("objective on the tape:", format(fit$opt$objective, digits = 12),
    " mcse:", format(fit$importance$mcse, digits = 4), " ess_min:",
    format(fit$importance$ess_min, digits = 5), "\n")
tpl <- fit$estimates
par <- fit$opt$par
pl_args <- length(formals(ns$imp_par_list))
pl <- if (pl_args >= 3L) {
  ns$imp_par_list(tpl, par, fit$frame[["map"]])
} else {
  ns$imp_par_list(tpl, par)
}
cat("betad the numeric path evaluates:", format(pl$betad, digits = 6),
    "\n  fit's betad:", format(tpl$betad, digits = 6),
    " log(5) =", format(log(5), digits = 6), "\n")
prop <- ns$imp_frozen_proposal(fit)
a <- as.matrix(prop$io$amat(pl))
lme <- function(v) max(v) + log(mean(exp(v - max(v))))
num <- -sum(apply(a, 1L, lme)) - sum(1 / 2 * log(2 * pi) + prop$plan$logdet_l)
cat("objective, numeric path:", format(num, digits = 12), "\n")
cat("difference:", format(num - fit$opt$objective, digits = 4), "\n")
