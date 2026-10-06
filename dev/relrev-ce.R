# Reviewer, item 3: conditional_effects(method = "predict") estimate__
# is the median of its simulated responses, on every family path. Per
# model: the display runs, lower <= estimate <= upper, whole counts for
# a discrete family (ndraws odd), and the estimate against an
# independent median: simulate(newdata = the display's grid) at the
# point estimate, 20001 draws per row. Seed 31 (data), 7 (display).
#   Rscript dev/relrev-ce.R <r5|r6> > dev/relrev-log/ce-<arm>.txt 2>&1
arm <- commandArgs(TRUE)[1]
lib <- if (arm == "r6") "C:/Users/adf44/source/r/rellib-r6" else
  "C:/Users/adf44/source/r/rellib-r5"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
set.seed(31)
n <- 300
d <- data.frame(x = runif(n, -1, 1), z = rnorm(n),
                g = factor(sample(letters[1:5], n, TRUE)))
d$yg <- 1 + d$x + rnorm(n)
d$yln <- exp(0.5 + d$x + rnorm(n, 0, 0.6))
d$yp <- rpois(n, exp(1 + d$x))
d$ynb <- rnbinom(n, mu = exp(1 + d$x), size = 2)
d$nt <- sample(5:15, n, TRUE)
d$yb <- rbinom(n, d$nt, plogis(d$x))
d$ybern <- rbinom(n, 1, plogis(d$x))
d$ytr <- 0.5 + abs(rnorm(n, d$x, 1))
d$cens <- ifelse(d$yg > 2, "right", "none"); d$ygc <- pmin(d$yg, 2)
d$yhu <- ifelse(runif(n) < 0.7, 0L, 1L + rpois(n, exp(0.5 + d$x)))
d$yzi <- ifelse(runif(n) < 0.6, 0L, rpois(n, exp(1 + d$x)))
d$ymix <- ifelse(runif(n) < 0.5, rnorm(n, -2 + d$x), rnorm(n, 2 + d$x))
d$ymi <- d$z; d$ymi[sample(n, 40)] <- NA
d$yo <- cut(d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
d$ycat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$yw <- rweibull(n, 2, exp(0.5 + 0.3 * d$x))
models <- list(
  gaussian = quote(frm(yg ~ x, family = gaussian(), data = d)),
  lognormal = quote(frm(yln ~ x, family = lognormal(), data = d)),
  weibull = quote(frm(yw ~ x, family = weibull(), data = d)),
  poisson = quote(frm(yp ~ x, family = poisson(), data = d)),
  negbinomial = quote(frm(ynb ~ x, family = negbinomial(), data = d)),
  binomial_trials = quote(frm(yb | trials(nt) ~ x, family = binomial(), data = d)),
  bernoulli = quote(frm(ybern ~ x, family = bernoulli(), data = d)),
  trunc_gaussian = quote(frm(ytr | trunc(lb = 0.5) ~ x, family = gaussian(), data = d)),
  cens_gaussian = quote(frm(ygc | cens(cens) ~ x, family = gaussian(), data = d)),
  hurdle_poisson = quote(frm(bf(yhu ~ x, hu ~ 1), family = hurdle_poisson(), data = d)),
  zi_poisson = quote(frm(bf(yzi ~ x, zi ~ 1), family = zero_inflated_poisson(), data = d)),
  mixture_gauss = quote(frm(bf(ymix ~ x), family = mixture(gaussian(), gaussian()), data = d)),
  mv_gauss_pois = quote(frm(mvbf(bf(yg ~ x), bf(yp ~ x)), data = d,
                            family = list(gaussian(), poisson()))),
  mi_model = quote(frm(bf(yg ~ x + mi(ymi)) + bf(ymi | mi() ~ z) + set_rescor(FALSE), data = d)),
  re_poisson = quote(frm(yp ~ x + (1 | g), family = poisson(), data = d)),
  cumulative = quote(frm(yo ~ x, family = cumulative(), data = d)),
  categorical = quote(frm(ycat ~ x, family = categorical(), data = d))
)
for (nm in names(models)) {
  fit <- tryCatch(suppressWarnings(eval(models[[nm]])),
                  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(fit, "err")) { cat(sprintf("%-16s FIT ERR %s\n", nm, substr(fit, 1, 100))); next }
  resp <- switch(nm, mv_gauss_pois = "yp", mi_model = "yg", NULL)
  set.seed(7)
  ce <- tryCatch(suppressWarnings(suppressMessages(
    conditional_effects(fit, effects = "x", method = "predict", ndraws = 1001,
                        resp = resp, resolution = 5))),
    error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(ce, "err")) { cat(sprintf("%-16s CE  ERR %s\n", nm, substr(gsub("\n", " ", ce), 1, 130))); next }
  df <- ce[[1]]
  est <- df$estimate__
  inside <- all(df$lower__ <= est & est <= df$upper__)
  whole <- all(est == round(est))
  # an independent median at the same grid rows, at the point estimate
  ind <- if (!is.null(resp)) structure("multivariate: skipped", class = "err") else tryCatch({
    nd <- df
    set.seed(8)
    s <- suppressWarnings(simulate(fit, nsim = 20001, newdata = nd,
                                   re_formula = NA))
    sm <- as.matrix(s)
    if (nrow(sm) != nrow(nd)) sm <- t(sm)
    apply(sm, 1, stats::median)
  }, error = function(e) structure(conditionMessage(e), class = "err"))
  cmp <- if (inherits(ind, "err")) paste("indep ERR", substr(ind, 1, 80)) else {
    sprintf("max|est - indep median| %.4g (est sd %.3g)", max(abs(est - ind)),
            mean(df$se__))
  }
  cat(sprintf("%-16s OK  est %s | inside %s whole %s | %s\n", nm,
              paste(format(est, digits = 4), collapse = " "), inside, whole, cmp))
}
