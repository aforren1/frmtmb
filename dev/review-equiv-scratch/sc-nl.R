# Short check of the nl sweep: 12 families on my own data (seed 5101,
# n = 300), Gamma and inverse.gaussian under the log link AND their
# default links, asym_laplace over three seeds. Outcome rule as the
# lane's: nl logLik equal to the linear fit's (1e-4 relative) honors the
# body; equal to the intercept-only fit drops it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ll <- function(expr) tryCatch(as.numeric(logLik(suppressWarnings(
  suppressMessages(expr)))), error = function(e)
    structure(NA_real_, err = substr(conditionMessage(e), 1, 110)))
one <- function(lab, fam, y, nl_rhs = "b0 + b1 * x", pars = "b0 + b1 ~ 1",
                lhs = "y", seed = 5101) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
  d$y <- y(d$x)
  a <- ll(frm(as.formula(paste(lhs, "~ x")), family = fam, data = d))
  b <- ll(frm(bf(as.formula(paste(lhs, "~", nl_rhs)), as.formula(pars),
                 nl = TRUE), family = fam, data = d))
  c0 <- ll(frm(as.formula(paste(lhs, "~ 1")), family = fam, data = d))
  tol <- function(v) 1e-4 * max(1, abs(v))
  out <- if (is.na(b)) paste("nl fails:", attr(b, "err")) else
    if (is.na(a)) paste("linear fails:", attr(a, "err")) else
    if (abs(a - b) < tol(a)) "honored" else
    if (abs(c0 - b) < tol(c0)) "DROPPED" else "differs"
  cat(sprintf("%-34s lin %10.4f  nl %10.4f  int %10.4f  %s\n", lab, a, b,
              c0, out))
}
pos <- function(x) exp(0.4 + 0.3 * x + rnorm(length(x), 0, 0.4))
one("Gamma(log)", Gamma(link = "log"), pos)
one("Gamma(inverse, default)", Gamma(), pos)
one("Gamma(inverse), body 1/(b0+b1*x)", Gamma(), pos,
    nl_rhs = "b0 + b1 * x")
one("inverse.gaussian(log)", inverse.gaussian(link = "log"), pos)
one("inverse.gaussian(1/mu^2, default)", inverse.gaussian(), pos)
for (s in c(5101, 5102, 5103)) {
  one(paste("asym_laplace seed", s), asym_laplace(),
      function(x) 0.3 + 0.5 * x + rnorm(length(x)), seed = s)
}
one("poisson", poisson(), function(x) rpois(length(x), exp(0.3 + 0.5 * x)))
one("binomial", binomial(), function(x) rbinom(length(x), 10,
    plogis(0.3 + 0.5 * x)), lhs = "y | trials(ntr)")
one("student", student(), function(x) 0.3 + 0.5 * x + rt(length(x), 4))
one("weibull", weibull(), pos)
one("zero_inflated_negbinomial", zero_inflated_negbinomial(), function(x)
  ifelse(runif(length(x)) < 0.2, 0L, rpois(length(x), exp(0.3 + 0.5 * x))))
one("hurdle_lognormal", hurdle_lognormal(), function(x)
  ifelse(runif(length(x)) < 0.2, 0, pos(x)))
one("von_mises", von_mises(), function(x)
  ((0.3 + 0.5 * x + rnorm(length(x), 0, 0.5) + pi) %% (2 * pi)) - pi)
ordy <- function(x) ordered(cut(0.3 + 0.5 * x + rlogis(length(x)),
                                c(-Inf, -0.5, 0.5, 1.5, Inf), labels = 1:4))
one("cumulative", cumulative(), ordy, nl_rhs = "b1 * x", pars = "b1 ~ 1")
one("hurdle_cumulative", hurdle_cumulative(), function(x) {
  o <- as.integer(ordy(x)); ordered(ifelse(runif(length(x)) < 0.2, 0L, o),
                                    levels = 0:4) },
  nl_rhs = "b1 * x", pars = "b1 ~ 1")
cat("compat now: Gamma", as.data.frame(frm_compat("nl", "Gamma"))$status,
    "| inverse.gaussian",
    as.data.frame(frm_compat("nl", "inverse.gaussian"))$status, "\n")
