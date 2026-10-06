# Reviewer: on test-se-check.R's gr(g, by = f) fixture (seed 11), where
# theta_2 and theta_3 are lost on one negative eigen-direction, which
# rows of the repaired covariance differ from the covariance with those
# two parameters held (their rows of Q removed)?
.libPaths(c("C:/Users/adf44/source/r/cifixrev-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] +
                       u[dd$g, 2] * dd$x, 1)
fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | gr(g, by = f))),
                            family = gaussian(), data = dd))
cat("theta:", format(signif(fit$estimates$theta, 4)), "\n")
print(ns$sdr_of(fit)$se_null)
Q <- ns$joint_precision(fit)
jc <- ns$get_joint_cov(fit)
keep <- setdiff(seq_len(nrow(Q)), jc$lost_pos)
Vind <- as.matrix(Matrix::solve(Q[keep, keep]))
r <- diag(jc$V)[keep] / diag(Vind)
lab <- ns$joint_coef_labels(fit, jc)[keep]
o <- order(-abs(log(r)))
print(data.frame(label = lab[o], repaired = diag(jc$V)[keep][o],
                 held = diag(Vind)[o], ratio = r[o])[1:8, ])
nd <- dd[c(1, 2, 9, 10), ]
lb <- frm_lp_basis(fit, newdata = nd)
A <- as.matrix(lb$A)
p <- lb$coef_pos
q_rep <- rowSums((A %*% jc$V[p, p]) * A)
Vh <- matrix(0, nrow(Q), nrow(Q))
Vh[keep, keep] <- Vind
q_held <- rowSums((A %*% Vh[p, p]) * A)
cat("prediction se, repaired:", format(signif(sqrt(q_rep), 5)), "\n")
cat("prediction se, held:    ", format(signif(sqrt(q_held), 5)), "\n")
