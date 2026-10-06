# Punch round 1, m2: what the exact AD Hessian costs against
# optimHess() on the reviewer's large models (dev/nanse-rev-cost.R's
# designs), and the fit's own time. Minimum of 3 timed repeats each.
#   Rscript dev/nanse-p1-hecost.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
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
d4 <- mk(3000, 50)
d5 <- mk(3000, 20, 30)
ctl <- frmtmb_control(check_se = "ignore")
models <- list(
  glm_f50 = function() frm(y ~ x + f, data = d4, family = poisson(),
                           control = ctl),
  glm_f300 = function() frm(y ~ x + f, data = d1, family = poisson(),
                            control = ctl),
  glmm_f20 = function() frm(y ~ x + f + (1 | g), data = d5,
                            family = poisson(), control = ctl),
  glmm_f100 = function() frm(y ~ x + f + (1 | g), data = d2,
                             family = poisson(), control = ctl),
  lmm_f400 = function() frm(yg ~ x + f + (1 | g), data = d3, control = ctl)
)
tm <- function(expr) {
  best <- Inf
  for (r in 1:3) {
    t0 <- proc.time()[[3]]
    force(eval.parent(substitute(expr)))
    best <- min(best, proc.time()[[3]] - t0)
  }
  best
}
for (nm in names(models)) {
  t_fit <- tm(f <- models[[nm]]())
  p <- f$opt$par
  t_oh <- tm(optimHess(p, f$obj$fn, f$obj$gr))
  t_he <- if (!length(f$obj$env$random)) tm(f$obj$he(p)) else NA
  t_gr <- tm(for (i in 1:20) f$obj$gr(p)) / 20
  cat(sprintf(paste0("%-10s p=%3d fit %.3fs optimHess %.3fs (%.2f of fit) ",
                     "he %s  one gradient %.4fs  nlminb gradients %s\n"),
              nm, length(p), t_fit, t_oh, t_oh / t_fit,
              if (is.na(t_he)) "n/a" else sprintf("%.3fs (%.2f)", t_he,
                                                  t_he / t_fit),
              t_gr, paste(f$opt$evaluations, collapse = "/")))
}
