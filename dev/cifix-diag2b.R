# Where does the base build's joint covariance go indefinite on the
# test-difference.R:490 fit, and which rows does the grid design load?
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r6"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(8)
n <- 120
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
fit <- suppressWarnings(
  frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),
      family = stats::gaussian(), data = d))
jc <- frmtmb:::get_joint_cov(fit)
print(table(jc$names))
cat("diag V by component:\n")
print(tapply(diag(jc$V), jc$names, function(v) format(range(v), digits = 4)))
Q <- frmtmb:::joint_precision(fit)
print(dim(Q)); print(table(rownames(Q)))
eq <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
cat("Q eigen range:", format(range(eq), digits = 6),
    " cond:", format(max(eq) / min(abs(eq)), digits = 4), "\n")
print(fit$frame$re_blocks |> lapply(function(b) c(b$term_label, b$dim,
                                                  b$n_levels)))
gx <- d$x[-1] - diff(d$x) / 2
A <- data.frame(x = gx, fac = factor("A", levels = levels(d$fac)))
lb <- frm_lp_basis(fit, newdata = A, extra_cov = TRUE)
cat("coef_pos components:\n"); print(table(jc$names[lb$coef_pos]))
cat("labels:\n"); print(lb$coef_names)
evs <- eigen(lb$V, symmetric = TRUE, only.values = TRUE)$values
cat("V[coef_pos] eigen:", format(evs, digits = 4), "\n")
cat("max |A| per column:\n")
print(signif(apply(abs(as.matrix(lb$A)), 2, max), 3))
cat("sdreport method:", fit$control$reml, "\n")
cat("se_lost:\n"); print(frmtmb:::sdr_of(fit)$se_lost)
