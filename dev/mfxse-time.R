# The time of marginaleffects' calls on a poisson (1 + x | g) fit, at 30 and
# 500 levels of 40 rows, one arm per process. Minimum of 3 rounds; each round
# also times a fixed control (20 solves of a 300 x 300 matrix), which reads the
# machine's load and must not move between the arms.
#   MFX_ARM=before|after Rscript dev/mfxse-time.R
# Data seed 7411.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
el <- function(e) {
  t0 <- proc.time()[["elapsed"]]
  force(e)
  proc.time()[["elapsed"]] - t0
}
set.seed(1)
A <- crossprod(matrix(rnorm(300 * 300), 300))
control <- function() for (i in 1:20) solve(A)
for (J in c(30L, 500L)) {
  set.seed(7411)
  m <- 40L
  n <- J * m
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
  d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
  fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
  invisible(quiet(get_vcov(fit)))
  nd <- d[seq(1, n, length.out = 12), ]
  tt <- replicate(3L, c(
    control = el(control()),
    predictions12 = el(quiet(predictions(fit, newdata = nd))),
    avg_predictions = el(quiet(avg_predictions(fit))),
    avg_slopes = el(quiet(avg_slopes(fit, variables = "x"))),
    avg_slopes_novcov = el(quiet(avg_slopes(fit, variables = "x",
                                            vcov = FALSE)))))
  cat(sprintf("TIME arm=%s J=%d rows=%d coefs=%d %s\n", mfx_arm, J, n,
              length(get_coef(fit)),
              paste(sprintf("%s=%.2f", rownames(tt), apply(tt, 1, min)),
                    collapse = " ")))
}
