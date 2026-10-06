# Punch round 1: tier 3's inputs on the gr(g, by = f) fixture of
# dev/nanse-rev-grby.R: the finite-difference noise, the spectrum, the
# thresholds, and what se_tier3() decides.
#   Rscript dev/nanse-p1-probe.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | gr(g, by = f))),
                            family = gaussian(), data = dd))
nm <- ns$outer_par_names(fit)
h <- ns$fit_outer_hessian(fit)
H <- h$H
E <- h$E
cat("diag H:", paste(signif(diag(H), 3), collapse = " "), "\n")
cat("row max E:", paste(signif(apply(E, 1, max), 3), collapse = " "), "\n")
D <- sqrt(abs(diag(H)))
Es <- E / outer(D, D)
cat("scaled noise row max:", paste(signif(apply(Es, 1, max), 3),
                                  collapse = " "), "\n")
cat("enorm", signif(sqrt(sum(Es^2)), 3), "\n")
e <- eigen(H / outer(D, D), TRUE)
cat("eig:", paste(signif(e$values, 3), collapse = " "), "\n")
an <- ns$cov_from_hessian(fit, H, E)
cat("tier", an$tier, "lost:", paste(names(an$lost), an$lost, sep = "=",
                                    collapse = " "), "\n")
cat("SE:", paste(sprintf("%s=%.4g", nm, suppressWarnings(sqrt(diag(an$shown)))),
                 collapse = " "), "\n")
