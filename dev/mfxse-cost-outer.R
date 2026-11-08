# The time of a sandwich carried into 1000 group-level effects (review m1):
# predictions() on 2 rows with vcov = ~ g, poisson (1 + x | g), 500 levels of
# 40 rows, data seed 7411 (dev/review-mfxse-cost.R's design). Also the inverse
# alone, dense Cholesky against Matrix::solve() of the sparse block.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
ns <- asNamespace("frmtmb")
set.seed(7411)
J <- 500L; m <- 40L; n <- J * m
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
nd2 <- d[c(1, n), ]
el <- function(e) { t0 <- proc.time()[["elapsed"]]; force(e); proc.time()[["elapsed"]] - t0 }
Q <- ns$joint_precision(fit)
r <- fit$obj$env$random
cat(sprintf("COST inner rows %d\n", length(r)))
cat(sprintf("COST chol2inv dense %.2f s\n",
            el(chol2inv(chol(as.matrix(Q[r, r]))))))
cat(sprintf("COST Matrix::solve sparse %.2f s\n",
            el(as.matrix(Matrix::solve(Q[r, r])))))
p <- NULL
cat(sprintf("COST predictions(2 rows, vcov = ~ g) %.2f s\n",
            el(p <- quiet(predictions(fit, newdata = nd2, vcov = ~ g)))))
p0 <- quiet(predictions(fit, newdata = nd2))
cat("COST se ~g:", sprintf("%.5f", p$std.error), " default:",
    sprintf("%.5f", p0$std.error), "\n")
cat(sprintf("COST vcov_cluster(~ g, full) %.2f s\n",
            el(vcov_cluster(fit, ~ g, type = "CR1", full = TRUE))))
cat(sprintf("COST get_vcov(vcov = ~ g) %.2f s\n",
            el(quiet(get_vcov(fit, vcov = ~ g)))))
cat(sprintf("COST get_vcov() default %.2f s\n", el(quiet(get_vcov(fit)))))
cat(sprintf("COST predictions(2 rows) default %.2f s\n",
            el(quiet(predictions(fit, newdata = nd2)))))
