# Reviewer probe for the families instrument: does frm() actually FIT
# each family the instrument counts (it checks the name only)? One
# intercept-plus-slope fit per family on simulated data, seed 20261013,
# n = 300, family passed by brms's name through frmtmb::brmsfamily().
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20261013)
n <- 300
x <- rnorm(n)
pos <- function() rgamma(n, 4, 4 / exp(0.5 + 0.3 * x))
cnt <- function() rpois(n, exp(0.7 + 0.3 * x))
p01 <- function() plogis(0.2 + 0.5 * x)
ord <- function() cut(0.5 * x + rlogis(n), c(-Inf, -1, 0, 1, Inf),
                      labels = FALSE)
zi <- function(v, p = 0.2) ifelse(runif(n) < p, 0, v)
bet <- function() rbeta(n, 2 * p01() * 5, 2 * (1 - p01()) * 5)
gen <- list(
  acat = ord, cratio = ord, sratio = ord, cumulative = ord,
  categorical = function() sample(1:3, n, TRUE),
  asym_laplace = function() 1 + x + rnorm(n),
  bernoulli = function() rbinom(n, 1, p01()),
  Beta = bet,
  binomial = function() rbinom(n, 10, p01()),
  exgaussian = function() 1 + x + rnorm(n) + rexp(n, 2),
  exponential = function() rexp(n, 1 / exp(0.5 + 0.3 * x)),
  Gamma = pos, gaussian = function() 1 + x + rnorm(n),
  geometric = cnt, hurdle_gamma = function() zi(pos()),
  hurdle_lognormal = function() zi(exp(rnorm(n, 0.3 * x))),
  hurdle_negbinomial = function() zi(pmax(1, cnt())),
  hurdle_poisson = function() zi(pmax(1, cnt())),
  inverse.gaussian = pos, lognormal = function() exp(rnorm(n, 0.3 * x)),
  negbinomial = cnt, poisson = cnt,
  shifted_lognormal = function() 0.2 + exp(rnorm(n, 0.3 * x)),
  skew_normal = function() 1 + x + rnorm(n),
  student = function() 1 + x + rt(n, 5),
  von_mises = function() {
    v <- 0.3 * x + rnorm(n, 0, 0.5); atan2(sin(v), cos(v)) },
  weibull = function() rweibull(n, 2, exp(0.5 + 0.3 * x)),
  zero_inflated_beta = function() zi(bet()),
  zero_inflated_binomial = function() zi(rbinom(n, 10, p01())),
  zero_inflated_negbinomial = function() zi(cnt()),
  zero_inflated_poisson = function() zi(cnt()),
  zero_one_inflated_beta = function() {
    v <- bet(); u <- runif(n); v[u < 0.1] <- 0; v[u > 0.9] <- 1; v },
  cox = function() rexp(n, exp(0.3 * x)),
  multinomial = NULL)
res <- character()
for (nm in names(gen)) {
  if (is.null(gen[[nm]])) next
  d <- data.frame(x = x, y = gen[[nm]](), tr = 10L)
  fo <- if (nm %in% c("binomial", "zero_inflated_binomial")) {
    bf(y | trials(tr) ~ x)
  } else bf(y ~ x)
  r <- tryCatch({
    f <- suppressWarnings(frm(fo, family = frmtmb::brmsfamily(nm),
                              data = d))
    sprintf("fit, code %d, logLik %.3f", f$opt$convergence,
            as.numeric(logLik(f)))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 120)))
  cat(sprintf("%-26s %s\n", nm, r))
}
# multinomial: brms's spelling, a matrix response with trials
Y <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, c(0.2, 0.3, 0.5))))
d <- data.frame(x = x, tr = 10L); d$Y <- Y
r <- tryCatch({
  f <- suppressWarnings(frm(bf(Y | trials(tr) ~ x),
                            family = frmtmb::brmsfamily("multinomial"),
                            data = d))
  sprintf("fit, code %d", f$opt$convergence)
}, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 160)))
cat(sprintf("%-26s %s\n", "multinomial", r))
