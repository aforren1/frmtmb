# Can frm_allfit() or a start show bernoulli's second maximum? The help sentence
# depends on it. Data as in dev/nlfix-bern.R.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
invisible(rlogis(n))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                   c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
set.seed(9205)
d$yg <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
d$yc <- rpois(n, 1.5 * exp(0.6 * d$x))
d$yzi <- ifelse(runif(n) < 0.25, 0L, rpois(n, 1.5 * exp(0.6 * d$x)))
d$yhu <- ifelse(runif(n) < 0.25, 0L, 1L + rpois(n, 1.5 * exp(0.6 * d$x)))
d$ypo <- rgamma(n, 3, 3 / (1.5 * exp(0.6 * d$x)))
d$ybi <- rbinom(n, 1, plogis(1.5 * exp(0.6 * d$x) - 2))
f <- bf(ybi ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
a <- suppressWarnings(frm(f, family = bernoulli(), data = d))
cat("default", as.numeric(logLik(a)), "conv", a$opt$convergence, "beta",
    a$estimates$beta, "\n")
for (sfm in c(TRUE, FALSE)) {
  af <- tryCatch(suppressWarnings(suppressMessages(
    frm_allfit(a, start_from_mle = sfm))), error = function(e) e)
  cat("\nfrm_allfit start_from_mle =", sfm, "\n")
  if (inherits(af, "error")) { cat("ERROR", conditionMessage(af), "\n"); next }
  print(af)
}
r <- suppressWarnings(frm(f, family = bernoulli(), data = d, start = list(
  beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.5))))
cat("\nstart b1 = 1:", as.numeric(logLik(r)), "conv", r$opt$convergence, "beta",
    r$estimates$beta, "\n")
