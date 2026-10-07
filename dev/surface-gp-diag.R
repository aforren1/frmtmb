# Lane surface, item 9: why frm_sample(fit) on an exact y ~ gp(x) fit
# does not move. Rebuilds the fit route's objective and init the way
# frm_sample() does and evaluates the joint density and its gradient
# there, beside the formula route's random init.
#
#   Rscript dev/surface-gp-diag.R lane|base > dev/surface-out/gp-diag-<arm>.txt
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
cat("outer estimates:\n")
print(fit$opt$par)
cat("template components and lengths:\n")
print(lengths(fit$frame$par_template))
cat("random:", fit$obj$env$random[1:3], "... of", length(fit$obj$env$random),
    "\n")
rp <- ns$sample_resolve_priors(fit, NULL, base = fit$prior, defaults = TRUE)
print(rp$effective)
obj <- ns$prior_augmented_obj(fit, rp$ri$entries)
lpb <- obj$env$last.par.best
cat("length(lpb)", length(lpb), "finite", sum(is.finite(lpb)), "\n")
f0 <- obj$env$f(lpb, order = 0)
g0 <- obj$env$f(lpb, order = 1)
cat("joint nll at the mode init:", format(f0, digits = 10), "\n")
cat("gradient at the init: finite", sum(is.finite(g0)), "of", length(g0),
    "| max abs", format(max(abs(g0[is.finite(g0)]))), "\n")
cat("names of non-finite gradient entries:",
    names(lpb)[!is.finite(g0)], "\n")
# a small step from the init along the gradient, as the first leapfrog
for (h in c(1e-8, 1e-6, 1e-4, 1e-2)) {
  x1 <- lpb - h * as.numeric(g0)
  cat(sprintf("step %.0e: nll %s\n", h,
              format(obj$env$f(x1, order = 0), digits = 10)))
}
# the field coordinates and their covariance at the mode
th <- fit$estimates$theta
cat("theta at the mode:", format(th, digits = 6), "\n")
b <- fit$estimates$b
cat("field: length", length(b), "range", format(range(b), digits = 4), "\n")
# the joint density's Hessian in the field block: its conditioning is
# what the leapfrog step size meets
H <- tryCatch(optimHess(lpb, function(p) obj$env$f(p, order = 0),
                        function(p) as.numeric(obj$env$f(p, order = 1))),
              error = function(e) NULL)
if (!is.null(H)) {
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  cat("joint Hessian eigenvalues: max", format(max(ev), digits = 4),
      "min", format(min(ev), digits = 4), "\n")
}
