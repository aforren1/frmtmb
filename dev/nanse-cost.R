# What a fit-time standard-error check would cost: per model, the fit
# time, the outer Hessian alone (optimHess on obj$fn/obj$gr, which is
# what sdreport computes first) and the full sdreport. Minimum of 3
# rounds each, interleaved.
#   Rscript dev/nanse-cost.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
set.seed(1)
sleep <- lme4::sleepstudy
epi <- brms::epilepsy
n <- 2000
dg <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:50, n, TRUE)))
dg$y <- rpois(n, exp(0.2 + 0.3 * dg$x + rnorm(50, 0, 0.5)[dg$g]))
dl <- data.frame(x = runif(300, 0, 3))
dl$y <- 2 * exp(-0.3 * dl$x) + rnorm(300, 0, 0.2)
models <- list(
  lm = function() frm(Reaction ~ Days, data = sleep),
  lmm_slope = function() frm(Reaction ~ Days + (Days | Subject), data = sleep),
  pois_glmm = function() frm(count ~ zAge + zBase * Trt + (1 | patient),
                             data = epi, family = poisson()),
  pois_glmm_2000 = function() frm(y ~ x + z + (1 | g), data = dg,
                                  family = poisson()),
  dist_sigma = function() frm(bf(Reaction ~ Days, sigma ~ Days), data = sleep),
  nl = function() frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                      data = dl, start = list(beta = c(2, -0.3))),
  negbin = function() frm(count ~ zAge + zBase * Trt + (1 | patient),
                          data = epi, family = negbinomial())
)
tm <- function(expr) {
  t0 <- proc.time()[[3]]
  force(expr)
  proc.time()[[3]] - t0
}
for (nm in names(models)) {
  tf <- th <- ts <- numeric(3)
  for (r in 1:3) {
    reps <- 0L
    t <- 0
    while (t < 1.2) {
      t <- t + tm(f <- suppressWarnings(models[[nm]]()))
      reps <- reps + 1L
    }
    tf[r] <- t / reps
    reps <- 0L
    t <- 0
    while (t < 1.2) {
      t <- t + tm(optimHess(f$opt$par, f$obj$fn, f$obj$gr))
      reps <- reps + 1L
    }
    th[r] <- t / reps
    reps <- 0L
    t <- 0
    while (t < 1.2) {
      t <- t + tm(frmtmb:::autoscale_sdreport(f))
      reps <- reps + 1L
    }
    ts[r] <- t / reps
  }
  cat(sprintf("%-15s p=%3d fit %.4fs hessian %.4fs (%.0f%%) sdreport %.4fs (%.0f%%)\n",
              nm, length(f$opt$par), min(tf), min(th),
              100 * min(th) / min(tf), min(ts), 100 * min(ts) / min(tf)))
}
