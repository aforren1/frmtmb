# Reviewer: the fit-time SE check's cost on models with many fixed
# effects, where optimHess() needs 2 gradients per parameter. Arms in
# one process, interleaved: check (default), ignore, and check again as
# the control (must read 1.0). Each arm grown past 1.2 s, minimum over
# rounds. Also: fit + summary() in both arms (the Hessian is said to be
# reused), and the size of the cached Hessian.
#   Rscript dev/nanse-rev-cost.R [rounds]
args <- commandArgs(trailingOnly = TRUE)
rounds <- if (length(args)) as.integer(args[1]) else 3L
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(11)
mk <- function(n, k, g = 0) {
  d <- data.frame(f = factor(sample(seq_len(k), n, TRUE)), x = rnorm(n))
  eta <- 0.3 + 0.2 * d$x + rnorm(k, 0, 0.3)[d$f]
  if (g) {
    d$g <- factor(sample(seq_len(g), n, TRUE))
    eta <- eta + rnorm(g, 0, 0.5)[d$g]
  }
  d$y <- rpois(n, exp(eta))
  d$yg <- eta + rnorm(n)
  d
}
d1 <- mk(6000, 300)
d2 <- mk(4000, 100, 40)
d3 <- mk(8000, 400, 60)
models <- list(
  glm_f300 = function(ctl) frm(y ~ x + f, data = d1, family = poisson(),
                               control = ctl),
  glmm_f100 = function(ctl) frm(y ~ x + f + (1 | g), data = d2,
                                family = poisson(), control = ctl),
  lmm_f400 = function(ctl) frm(yg ~ x + f + (1 | g), data = d3,
                               control = ctl)
)
on <- frmtmb_control()
off <- frmtmb_control(check_se = "ignore")
per_call <- function(f, ctl, then = NULL) {
  reps <- 0L
  t0 <- proc.time()[[3]]
  while (proc.time()[[3]] - t0 < 1.2) {
    ft <- suppressWarnings(f(ctl))
    if (!is.null(then)) invisible(suppressWarnings(then(ft)))
    reps <- reps + 1L
  }
  (proc.time()[[3]] - t0) / reps
}
summ <- function(ft) summary(ft)
for (nm in names(models)) {
  ft <- suppressWarnings(models[[nm]](on))
  h <- ft$cache$hessian_fixed
  cat(sprintf("%-10s outer pars %d  code %d  cached Hessian %s\n", nm,
              length(ft$opt$par), ft$opt$convergence,
              format(object.size(as.list(ft$cache)), units = "MB")))
  a <- b <- cc <- sa <- sb <- numeric(rounds)
  for (r in seq_len(rounds)) {
    a[r] <- per_call(models[[nm]], on)
    b[r] <- per_call(models[[nm]], off)
    cc[r] <- per_call(models[[nm]], on)
    sa[r] <- per_call(models[[nm]], on, summ)
    sb[r] <- per_call(models[[nm]], off, summ)
  }
  cat(sprintf(paste0("%-10s fit: check %.3fs ignore %.3fs ratio %.3f ",
                     "control %.3f | fit+summary: check %.3fs ignore ",
                     "%.3fs ratio %.3f\n"),
              nm, min(a), min(b), min(a) / min(b), min(cc) / min(a),
              min(sa), min(sb), min(sa) / min(sb)))
}
