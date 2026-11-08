# Probe: the cost of one get_predict() and of the marginaleffects calls on a
# mixed fit, with the clock checked (data seed 7411).
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(7411)
J <- as.integer(Sys.getenv("MFX_J", "30")); m <- 40; n <- J * m
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
tm <- function(lab, e, reps = 1L) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) r <- force(e())
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("%-28s %.4f s per call (%d calls, %.2f s)\n", lab, el / reps,
              reps, el))
  invisible(r)
}
b0 <- fit$estimates$b
tm("get_predict same fit", function() get_predict(fit, newdata = d), 200)
tm("get_predict moved b", function() {
  f <- fit
  f$estimates$b <- b0 + rnorm(length(b0), 0, 1e-6)
  get_predict(f, newdata = d)
}, 200)
tm("frm_linpred moved b", function() {
  f <- fit
  f$estimates$b <- b0 + rnorm(length(b0), 0, 1e-6)
  frm_linpred(f, newdata = d, type = "response")
}, 200)
r <- tm("avg_slopes", function() quiet(avg_slopes(fit, variables = "x")))
print(as.data.frame(r)[, c("term", "estimate", "std.error")])
r <- tm("avg_predictions", function() quiet(avg_predictions(fit)))
print(as.data.frame(r)[, c("estimate", "std.error")])
r <- tm("predictions 12 rows", function() {
  quiet(predictions(fit, newdata = d[seq(1, n, length.out = 12), ]))
})
