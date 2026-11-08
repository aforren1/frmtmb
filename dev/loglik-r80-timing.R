# What the split costs per logLik() call: logLik() now resolves the
# fit's prior and evaluates its log density at the estimates. Arms are
# interleaved in one process, each block grown past 1.2 s, minimum of 5
# rounds; the control is logLik() on the fit with no prior (no resolve),
# and logLik_penalized() on the no-prior fit is the arithmetic floor.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(4101)
n <- 200; ng <- 20
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = n / ng)))
d$y <- 1 + 0.5 * d$x + rnorm(ng)[d$g] + rnorm(n)
f0 <- frm(y ~ x + (1 | g), data = d)
pr <- set_prior("normal(0, 1)", class = "b") +
  set_prior("student_t(3, 0, 1)", class = "sd") +
  set_prior("normal(1, 2)", class = "Intercept")
fp <- frm(y ~ x + (1 | g), data = d, prior = pr)
arms <- list(
  control_noprior = function() logLik(f0),
  prior = function() logLik(fp),
  penalized_prior = function() logLik_penalized(fp))
reps <- vapply(arms, function(f) {
  k <- 1L
  repeat {
    t <- system.time(for (i in seq_len(k)) f())[["elapsed"]]
    if (t > 1.2) break
    k <- k * 2L
  }
  k
}, 1L)
res <- matrix(NA_real_, 5, length(arms), dimnames = list(NULL, names(arms)))
for (r in 1:5) {
  for (a in names(arms)) {
    f <- arms[[a]]
    k <- reps[[a]]
    res[r, a] <- system.time(for (i in seq_len(k)) f())[["elapsed"]] / k
  }
}
cat("reps per block:", paste(names(reps), reps, sep = "=", collapse = " "),
    "\n")
cat(sprintf("%-16s min %.1f us  median %.1f us\n", colnames(res),
            apply(res, 2, min) * 1e6, apply(res, 2, stats::median) * 1e6),
    sep = "")
# and against one fit, the unit the call is made once per
tf <- system.time(frm(y ~ x + (1 | g), data = d, prior = pr))[["elapsed"]]
cat(sprintf("one MAP fit: %.3f s\n", tf))
