# Punch round 1, m3: which warnings the frm_sample() calls of
# extensions/frmtmb.sample/tests/testthat/test-ordinal-mixture-draws.R
# raise (the fits and draws there, at the same seeds), so the test can
# name them in allow_warnings() rather than suppress every warning.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
omd_data <- function(seed, n = 200) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < 0.2, 0L, d$y)
  d
}
grab <- function(label, expr) {
  w <- character(0)
  withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  cat("==", label, length(w), "warnings\n")
  for (m in w) cat("   WARN", substr(gsub("\n", " ", m), 1, 220), "\n")
}
d <- omd_data(71)
st <- "student_t(3, 0, 2.5)"
nb <- set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
  set_prior("normal(0, 2)", class = "b", dpar = "mu2")
fits <- list(
  none = frm(bf(y ~ x, disc1 ~ 0 + z),
             family = mixture(cumulative(), sratio()), data = d,
             prior = nb +
               set_prior(st, class = "Intercept", dpar = "mu1") +
               set_prior(st, class = "Intercept", dpar = "mu2")),
  mu = frm(bf(y ~ x),
           family = mixture(cumulative(), acat(threshold = "sum_to_zero"),
                            order = "mu"),
           data = d, prior = nb + set_prior(st, class = "Intercept")),
  hgr = frm(bf(yh | thres(gr = g) ~ x), family = hurdle_cumulative(),
            data = d, prior = set_prior(st, class = "Intercept")),
  hcs = frm(bf(yh ~ cs(x)), family = hurdle_cumulative("probit"),
            data = d,
            prior = set_prior(st, class = "Intercept") +
              set_prior("normal(0, 1)", class = "b")))
for (k in names(fits)) {
  grab(k, frm_sample(fits[[k]], chains = 1, iter = 300, refresh = 0,
                     seed = 3))
}
set.seed(20261005)
d <- data.frame(z = rnorm(200))
u <- stats::rlogis(200) / exp(0.4 * d$z)
d$y <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
fit <- frm(bf(y ~ 1, disc ~ 0 + z), family = cumulative(), data = d,
           prior = set_prior("student_t(3, 0, 2.5)", class = "Intercept"))
grab("noloc", frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3))
