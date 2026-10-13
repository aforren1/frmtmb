# Reviewer of lane seflags, re-check: the cost of se_design_null()'s QR
# on a wide design just under se_design_flops (1e10). y ~ f + x + s(x)
# with a 500-level factor on 36000 rows (n p^2 = 9.4e9): the smooth's
# fixed column is aliased with x, so a direction is removed and the QR
# runs. Times se_check() on a fit made with the check off, the QR alone,
# and the check's evaluation counts.
#   Rscript dev/seflags-rev2-qrcost.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R")
set.seed(1)
n <- 36000
d <- data.frame(f = factor(rep(1:500, length.out = n)), x = rnorm(n))
d$y <- rnorm(500)[d$f] + 0.5 * d$x + sin(2 * d$x) + rnorm(n)
t0 <- system.time(fit <- suppressWarnings(suppressMessages(
  frm(y ~ f + x + s(x, k = 10), data = d,
      control = frmtmb_control(check_se = "ignore")))))[["elapsed"]]
cat("fit (check off)", t0, "s | outer", length(fit$opt$par), "| X",
    paste(dim(fit$frame$linpreds[[1]]$X), collapse = " x "), "\n")
n_fn <- 0L; n_gr <- 0L
fn0 <- fit$obj$fn; gr0 <- fit$obj$gr
fit$obj$fn <- function(...) { n_fn <<- n_fn + 1L; fn0(...) }
fit$obj$gr <- function(...) { n_gr <<- n_gr + 1L; gr0(...) }
tc <- system.time(invisible(suppressWarnings(suppressMessages(
  ns$se_check(fit, frmtmb_control(check_se = "warning"))))))[["elapsed"]]
cat("se_check alone", tc, "s | fn", n_fn, "gr", n_gr, "| lost",
    substr(lost_of(fit), 1, 120), "\n")
if (exists("se_design_null", ns)) {
  tq <- system.time(N <- ns$se_design_null(fit))[["elapsed"]]
  cat("se_design_null alone", tq, "s | columns", NCOL(N), "\n")
}
cat("peak memory (MB, gc max used):",
    sum(gc()[, ncol(gc()) - 1L]), "\n")
