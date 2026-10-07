# Lane surface, item 9, fourth step: the first warmup iterations from
# the mode, and the joint density far from it, where a large first step
# lands.
#
#   Rscript dev/surface-gp-diag4.R lane|base > dev/surface-out/gp-diag4-<arm>.txt
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "| frmtmb", as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb.sample")
set.seed(5)
n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
rp <- suppressMessages(ns$sample_resolve_priors(fit, NULL, base = fit$prior,
                                                defaults = TRUE))
obj <- ns$prior_augmented_obj(fit, rp$ri$entries)
lpb <- obj$env$last.par.best
nm <- names(lpb)
sf <- suppressWarnings(tmbstan::tmbstan(
  obj, chains = 1, iter = 300, refresh = 0, seed = 4,
  init = list(array(lpb, dim = length(lpb)))))
sp <- rstan::get_sampler_params(sf, inc_warmup = TRUE)[[1]]
print(round(sp[1:25, ], 4))
a <- rstan::extract(sf, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
cat("distinct draws over all 300 iterations:",
    nrow(unique(round(a[, 1:5], 12))), "\n")
cat("\njoint nll far from the mode, per coordinate group:\n")
for (k in unique(nm)) {
  for (dlt in c(-50, -20, -5, 5, 20, 50)) {
    x <- lpb
    x[nm == k] <- x[nm == k] + dlt
    f <- obj$env$f(x, order = 0)
    g <- obj$env$f(x, order = 1)
    cat(sprintf("  %-6s %+4d: nll %-14s gradient finite %d of %d\n", k, dlt,
                format(f, digits = 8), sum(is.finite(g)), length(g)))
  }
}

cat("\nthe point the chain moved to:\n")
u <- unique(round(a, 12))
x2 <- as.numeric(a[nrow(a), seq_along(nm)])
names(x2) <- nm
print(signif(x2[nm != "b"], 6))
cat("field range", format(range(x2[nm == "b"]), digits = 6), "\n")
cat("nll there:", obj$env$f(x2, order = 0), "\n")
# the joint's pieces at that point: likelihood alone, without the field
# prior, read from the plain objective
nll0 <- frmtmb::build_objective(fit$frame)
pp <- fit$estimates
pp$beta[] <- x2[nm == "beta"]; pp$betad[] <- x2[nm == "betad"]
pp$b[] <- x2[nm == "b"]; pp$theta[] <- x2[nm == "theta"]
cat("bare joint nll (no prior) there:", nll0(pp), "\n")
blk <- fit$frame$re_blocks[[1]]
gpn <- frmtmb:::covstruct_registry[["gp"]]$nll
cat("gp field log density there:",
    gpn(x2[nm == "b"], x2[nm == "theta"], blk), "\n")
th <- x2[nm == "theta"]
k <- length(x2[nm == "b"])
R <- frmtmb:::gp_corr(th, blk)
cat("log-scale form there:",
    mvtnorm_like <- sum(RTMB::dmvnorm(t(matrix(x2[nm == "b"] / exp(th[1]),
                                              nrow = 1)), 0, R,
                                     log = TRUE)) - k * th[1], "\n")
