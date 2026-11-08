# Item 3: test-gp-start.R's "y's and x's units" designs (seeds 4140, 1, 2,
# 20261102). At each scale: which route the outer Hessian took, the marginal
# against the z tape at the fit's point, and the fitted standard errors from
# each route's Hessian. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
orig <- get("gp_marginal_outer_gr", envir = ns)
use <- function(new) {
  unlockBinding("gp_marginal_outer_gr", ns)
  assign("gp_marginal_outer_gr", if (new) orig else function(fit, at) NULL,
         envir = ns)
  lockBinding("gp_marginal_outer_gr", ns)
}
designs <- list(
  list(seed = 4140, n = 40L, to = 6, f = function(x) sin(6 * x), s = 1e-3),
  list(seed = 1, n = 80L, to = 10, s = 1e-3, f = function(x) 2 * sin(3 * x) + x),
  list(seed = 2, n = 60L, to = 10, s = 0.1, f = function(x) sin(4 * x) + 0.5 * x),
  list(seed = 20261102, n = 60L, to = 6, f = function(x) sin(2 * x), s = 1e-3))
for (ds in designs) {
  set.seed(ds$seed)
  x <- round(stats::runif(ds$n, 0, ds$to), 3)
  y <- ds$f(x) + ds$s * stats::rnorm(ds$n)
  nd <- seq(min(x) + 0.05, max(x) - 0.07, length.out = 9L)
  for (sc in list(c(1, 1), c(1000, 1), c(1e-3, 1), c(1, 1000), c(1, 1e-3))) {
    se <- list()
    for (new in c(FALSE, TRUE)) {
      use(new)
      f <- suppressWarnings(frm(y ~ gp(x), data = data.frame(x = x * sc[2], y = y * sc[1])))
      ft <- suppressWarnings(fitted(f, newdata = data.frame(x = nd * sc[2])))
      se[[as.character(new)]] <- ft[, 2L] / sc[1]
    }
    use(TRUE)
    at <- f$opt$par
    s <- ns$gp_marginal_setup(f$frame, max_rows = Inf)
    p <- list(beta = f$estimates$beta[s$idx], theta = f$estimates$theta,
              betad = f$estimates$betad)
    v <- ns$gp_marginal_nll(s)(p)
    cat(sprintf("SD seed %d scale y %g x %g route %s valdiff %.2e | fitted se new/old max %.4f min %.4f | sigma %.3g units %s\n",
                ds$seed, sc[1], sc[2], is.function(orig(f, at)),
                abs(v - f$opt$objective) / max(1, abs(f$opt$objective)),
                max(se[["TRUE"]] / se[["FALSE"]]), min(se[["TRUE"]] / se[["FALSE"]]),
                exp(f$estimates$betad[[1]]) / sc[1],
                paste(signif(f$par_units %||% 1, 3), collapse = ",")))
  }
}
cat("DONE\n")
