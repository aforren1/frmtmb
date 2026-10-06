# Reviewer, final check: which Hessian and noise does the lane's analysis
# read on the RB4 fixture (seed 11)? cov_from_hessian() on
# fit_outer_hessian()'s H and E, against what sdr_of() reports.
#   Rscript dev/nanse-rev3-rb4-path.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(11)
n <- 480
d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                g2 = factor(rep(1:6, length.out = n)))
d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$f] + rnorm(n)
fit <- suppressWarnings(suppressMessages(frm(y ~ x + f + (1 + x | g2), data = d)))
cat("par_units all 1:", is.null(fit$par_units) || all(fit$par_units == 1), "\n")
h <- ns$fit_outer_hessian(fit)
an <- ns$cov_from_hessian(fit, h$H, h$E)
cat("cov_from_hessian(fit_outer_hessian) tier", an$tier, "lost", length(an$lost),
    paste(head(names(an$lost), 4), collapse = ","), "\n")
sdr <- suppressWarnings(ns$sdr_of(fit))
cat("sdr_of lost", length(sdr$se_lost), "\n")
hc <- fit$cache$hessian_fixed
cat("cached hessian_fixed:", !is.null(hc), " identical H:", identical(hc$H, h$H), "\n")
th <- which(ns$outer_par_names(fit) == "theta_3")
cat("theta_3 row max", max(abs(h$H[th, ])), " row noise max", max(h$E[th, ]),
    " diag", h$H[th, th], "\n")
