# Lane seflags, item 5: test-bcm-latent-mixtures.R's Malingering_2 fit,
# whose component-2 precision runs to the binomial limit. What does the
# SE check read: the spectrum, the corrected curvature of the smallest
# direction, the verdict, the warnings?
#   Rscript dev/seflags-bcmmal.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/tests/testthat"
env <- new.env(parent = asNamespace("frmtmb"))
ex <- parse(file.path(wt, "test-bcm-latent-mixtures.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
cat("lib", find.package("frmtmb"), "| R", R.home(), "\n")
d <- env$bcm_malingering_data()
w <- character()
fit <- withCallingHandlers(
  frm(bf(k | trials(n) ~ 1), family = mixture(beta_binomial, beta_binomial),
      data = d, start = env$bcm_mix_start(c(0.55, 0.95), phi = c(15, 15),
                                          theta = 0)),
  warning = function(x) { w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning") })
cat("code", fit$opt$convergence, "warnings", length(w), "\n")
for (x in w) cat("  W:", substr(x, 1, 300), "\n")
p <- fit$opt$par
print(p)
H <- fit$obj$he(p); H <- (H + t(H)) / 2
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("diag H:", signif(diag(H), 3), "\n")
cat("ev:", signif(e$values, 3), "\n")
cat("grad:", signif(fit$obj$gr(p), 3), "\n")
k <- length(e$values)
c0 <- if (exists("se_ridge_curvature", ns)) {
  ns$se_ridge_curvature(fit, p, seq_along(p), e, D, k, rep(TRUE, k))
} else NA
cat("corrected curvature of the smallest:", c0, "\n")
print(suppressWarnings(ns$sdr_of(fit))$se_lost)
