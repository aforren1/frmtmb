# Short check of 0.82.0 after the fix pass, probe 1: boundfam's
# consolidation fixes. (a) jc_nonest()'s move test on zi ~ x / hu ~ x over
# seeds and at new x; (c) jc_nan_lost() on ordinary fits with a parameter
# that has no standard error. usage: rel082-short-probe1.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
quiet_fit <- function(expr) {
  withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"),
                      message = function(m) invokeRestart("muffleMessage"))
}
fitted_report <- function(label, fit, newdata = NULL) {
  ws <- character()
  ft <- tryCatch(withCallingHandlers(
    if (is.null(newdata)) fitted(fit) else fitted(fit, newdata = newdata),
    warning = function(w) {
      ws <<- c(ws, substr(conditionMessage(w), 1, 70))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(ft, "error")) {
    cat(sprintf("%-34s ERROR %s\n", label, conditionMessage(ft)))
    return(invisible())
  }
  se <- ft[, "Est.Error"]
  cat(sprintf("%-34s rows %4d NaN %4d medSE %.4g warn %d %s\n", label,
              length(se), sum(!is.finite(se)), stats::median(se, na.rm = TRUE),
              length(ws), if (length(ws)) ws[1] else ""))
}
zib <- function(seed) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
  d
}
cat("\n## (a) zi ~ x, no excess zeros, seeds 1-12\n")
for (s in 1:12) {
  d <- zib(s)
  f <- quiet_fit(frm(bf(y | trials(N) ~ x, zi ~ x),
                     family = zero_inflated_binomial(), data = d))
  b <- fixef(f)[c("zi_Intercept", "zi_x"), "Estimate"]
  fitted_report(sprintf("zib s%-2d zi_I %7.2f zi_x %6.2f", s, b[1], b[2]), f)
}
cat("\n## (a) hu ~ x, no zeros, seeds 1-8\n")
for (s in 1:8) {
  set.seed(100 + s)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rpois(300, exp(1.5 + 0.3 * d$x))
  d$y[d$y == 0] <- 1L
  f <- quiet_fit(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = d))
  b <- fixef(f)[c("hu_Intercept", "hu_x"), "Estimate"]
  fitted_report(sprintf("hu  s%-2d hu_I %7.2f hu_x %6.2f", s, b[1], b[2]), f)
}
cat("\n## (a) zi ~ x seed 7 at new x: the undetermined slope reads there\n")
d <- zib(7)
f7 <- quiet_fit(frm(bf(y | trials(N) ~ x, zi ~ x),
                    family = zero_inflated_binomial(), data = d))
print(fixef(f7)[, 1:2])
for (x in c(0, 3, 10, 30, 100)) {
  fitted_report(sprintf("zib s7 newdata x = %g", x), f7,
                data.frame(x = x, N = 10))
}
cat("\n## (a) zi on a scaled x (x * 1e-4)\n")
d$xs <- d$x * 1e-4
fs <- quiet_fit(frm(bf(y | trials(N) ~ x, zi ~ xs),
                    family = zero_inflated_binomial(), data = d))
print(fixef(fs)[, 1:2])
fitted_report("zib s7 zi ~ xs", fs)
fitted_report("zib s7 zi ~ xs, newdata xs = 0.01", fs,
              data.frame(x = 0, xs = 0.01, N = 10))

cat("\n## (c) ordinary fits with a parameter that has no SE\n")
set.seed(5)
g <- factor(rep(1:10, each = 10))
x <- stats::rnorm(100)
dd <- data.frame(g = g, x = x, y = 1 + 0.5 * x + stats::rnorm(100))
fitted_report("gaussian (1|g) sd at 0",
              quiet_fit(frm(y ~ x + (1 | g), data = dd)))
fitted_report("gaussian (1|g) sd at 0, re_formula NA",
              quiet_fit(frm(y ~ x + (1 | g), data = dd)), NULL)
dd$yc <- 1 + 0.5 * x + (0.8 * stats::rnorm(10))[g] * (1 + x) +
  stats::rnorm(100, 0, 0.3)
fitted_report("gaussian (1+x|g), cor near 1",
              quiet_fit(frm(yc ~ x + (1 + x | g), data = dd)))
fitted_report("slope at a uniform prior bound",
              quiet_fit(frm(y ~ x, data = dd, prior = set_prior(
                "uniform(0, 0.2)", class = "b", coef = "x"))))
dd$k <- stats::rpois(100, exp(1 + 0.3 * x))
fitted_report("negbinomial on poisson data",
              quiet_fit(frm(k ~ x, data = dd, family = negbinomial())))
fitted_report("negbinomial on poisson, (1|g)",
              quiet_fit(frm(k ~ x + (1 | g), data = dd,
                            family = negbinomial())))
fitted_report("student on gaussian data",
              quiet_fit(frm(y ~ x, data = dd, family = student())))
dd$n <- 10L
dd$b <- stats::rbinom(100, 10, stats::plogis(0.2 * x))
fitted_report("beta_binomial on binomial data",
              quiet_fit(frm(b | trials(n) ~ x, data = dd,
                            family = beta_binomial())))
fitted_report("asym_laplace y ~ x",
              quiet_fit(frm(y ~ x, data = dd, family = asym_laplace())))
dd$sep <- as.integer(x > 0)
fitted_report("bernoulli separated",
              quiet_fit(frm(sep ~ x, data = dd, family = bernoulli())))
dd$z <- stats::rnorm(100)
fitted_report("bernoulli separated + z",
              quiet_fit(frm(sep ~ x + z, data = dd, family = bernoulli())))
