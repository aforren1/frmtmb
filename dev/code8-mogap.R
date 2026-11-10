# test-mo-coding.R's gap fit at cor = 1 (seed 23): conditions and the
# covariance get_vcov() returns, on each build.
source("dev/code8-common.R")
set.seed(23)
n <- 300
d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                g = factor(sample(1:10, n, TRUE)))
d$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d$g]) * d$x + 0.5 * d$z + rnorm(n)
r <- rfit(y ~ mo(x) + z + (mo(x) | g), data = d)
cat("code", r$f$opt$convergence, r$f$opt$message, "\n")
for (w in r$w) cat("WARN:", substr(w, 1, 300), "\n")
for (m in r$m) cat("MSG:", substr(m, 1, 300), "\n")
sdr <- quiet(fns$sdr_of(r$f))
cat("lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
print(data.frame(name = fns$outer_par_names(r$f), est = r$f$opt$par,
                 se = sqrt(diag(sdr$cov.fixed))), digits = 4)
V <- withCallingHandlers(marginaleffects::get_vcov(r$f), warning = function(w) {
  cat("GETVCOV WARN:", substr(conditionMessage(w), 1, 200), "\n")
  invokeRestart("muffleWarning") })
print(signif(V, 3))
