# lane tmbstan, probe 07: can rstan::grad_log_prob() catch a broken
# build directly, with no chain, no seed and no platform variance?
#
# The mechanism, from the reviewer. The two generated log_prob_impl
# overloads are selected by stan::require_not_st_var<VecR> (the
# PATCHED one) and stan::require_st_var<VecR> (the UNPATCHED one), so
# on a broken build the double-precision entry point reads the model
# and the reverse-mode entry point reads a standard normal.
# grad_log_prob() is reverse-mode. It should therefore return -u on a
# broken build, whatever the model is.
#
# This machine's tmbstan is CLEAN and no broken one may be installed,
# so the broken arm is not sampled. It is taken from the closed form
# the defect has: the standard normal kernel is -sum(u^2)/2 and its
# gradient is -u. That is not invented. It is the column
# dev/prior-dropping-investigation.md recorded from the affected
# container, where stan_lp came back as -mu^2/2 exactly and stan_gr as
# -mu.
#
# Sign convention: obj$fn is the NEGATIVE log posterior, Stan's
# log_prob is the log posterior, so the comparison is
# grad_log_prob(sf, u) against -obj$gr(u).
#
# SEED 4021 for the data, sampler seed 11.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-tmbstan/dev/stan-cache")
library(frmtmb)
library(frmtmb.sample)

set.seed(4021)
n <- 80L
dd <- data.frame(x = stats::rnorm(n))
dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)

cat("== does as_tmbstan() give a stanfit grad_log_prob can read? ==\n")
sf <- suppressWarnings(suppressMessages(
  as_tmbstan(fit, chains = 1L, iter = 20L, warmup = 10L, refresh = 0,
             seed = 11)))
cat("class:", class(sf), "\n")
np <- rstan::get_num_upars(sf)
cat("num upars:", np, "\n")
cat("fit obj par length:", length(fit$obj$par), "\n")

# points that are NOT the mode and not each other, so a gradient that
# ignores the data cannot coincide by luck
set.seed(4022)
pts <- list(rep(0, np), fit$obj$par + 0.3, fit$obj$par - 0.7,
            stats::rnorm(np))
names(pts) <- c("origin", "mode+0.3", "mode-0.7", "random")

cat("\n== clean arm: this installation ==\n")
rows <- list()
for (nm in names(pts)) {
  u <- as.numeric(pts[[nm]])
  g_stan <- rstan::grad_log_prob(sf, u)
  g_obj <- -as.numeric(fit$obj$gr(u))
  g_bad <- -u                      # the closed form of the defect
  scale <- max(abs(g_obj))
  rows[[nm]] <- data.frame(
    point = nm,
    max_abs_g_obj = scale,
    d_clean = max(abs(g_stan - g_obj)),
    d_clean_rel = max(abs(g_stan - g_obj)) / scale,
    d_broken_rel = max(abs(g_bad - g_obj)) / scale,
    identical = isTRUE(all.equal(g_stan, g_obj, tolerance = 0)),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
print(res, row.names = FALSE, digits = 6)

cat("\n== is the clean agreement exact, or only close? ==\n")
u <- as.numeric(fit$obj$par + 0.3)
g_stan <- rstan::grad_log_prob(sf, u)
g_obj <- -as.numeric(fit$obj$gr(u))
cat("identical():", identical(g_stan, g_obj), "\n")
cat("max ulp gap:",
    max(abs(g_stan - g_obj)) / max(.Machine$double.eps * abs(g_obj)), "\n")
cat("g_stan:", paste(format(g_stan, digits = 15), collapse = " "), "\n")
cat("g_obj :", paste(format(g_obj, digits = 15), collapse = " "), "\n")

cat("\n== the value entry point, for contrast ==\n")
# log_prob() is the double-precision entry, which on a BROKEN build
# reads the correct density. If so it cannot discriminate and the
# assertion must be on the gradient.
lp <- rstan::log_prob(sf, u)
cat("log_prob:", format(lp, digits = 12),
    " -obj$fn:", format(-as.numeric(fit$obj$fn(u)), digits = 12), "\n")

cat("\n== separation ==\n")
cat("clean relative gap, worst point :",
    format(max(res$d_clean_rel), digits = 4), "\n")
cat("broken relative gap, worst point:",
    format(max(res$d_broken_rel), digits = 4), "\n")
cat("broken relative gap, BEST point :",
    format(min(res$d_broken_rel), digits = 4), "\n")
