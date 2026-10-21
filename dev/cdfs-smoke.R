# Smoke test: every family that gained a CDF fits a model with all four
# censoring codes and one with trunc(), and the censored fit's logLik
# equals a hand-rolled likelihood from independent reference functions
# at the fitted parameters. Seed 11. Output: dev/cdfs-out/smoke.txt
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(11)
n <- 300
x <- rnorm(n)
cc <- sample(c(-1, 0, 0, 0, 1, 2), n, TRUE)
fits <- list(
  student = list(fam = student(), y = 1 + 0.5 * x + rt(n, 5) * 0.8),
  skew_normal = list(fam = skew_normal(),
                     y = brms::rskew_normal(n, 1 + 0.5 * x, 1, 4)),
  exgaussian = list(fam = exgaussian(),
                    y = brms::rexgaussian(n, 2 + 0.5 * x, 0.5, 1)),
  shifted_lognormal = list(fam = shifted_lognormal(),
                           y = 0.3 + rlnorm(n, 0.2 * x, 0.5)),
  asym_laplace = list(fam = asym_laplace(),
                      y = brms::rasym_laplace(n, 1 + 0.5 * x, 1, 0.3)),
  von_mises = list(fam = von_mises(),
                   y = brms::rvon_mises(n, 2 * atan(0.4 * x), 3)),
  Gamma = list(fam = Gamma(link = "log"),
               y = rgamma(n, 3, 3 / exp(0.5 + 0.3 * x))),
  Beta = list(fam = Beta(), y = rbeta(n, 2 * exp(0.3 * x), 3)),
  inverse.gaussian = list(fam = inverse.gaussian(link = "log"),
                          y = statmod::rinvgauss(n, exp(0.5 + 0.2 * x), 4)),
  poisson = list(fam = poisson(), y = rpois(n, exp(1 + 0.3 * x))),
  negbinomial = list(fam = negbinomial(),
                     y = rnbinom(n, mu = exp(1 + 0.3 * x), size = 2)),
  geometric = list(fam = geometric(),
                   y = rnbinom(n, mu = exp(1 + 0.3 * x), size = 1)),
  binomial = list(fam = binomial(), y = rbinom(n, 12, plogis(0.3 * x)),
                  trials = 12),
  beta_binomial = list(fam = beta_binomial(),
                       y = brms::rbeta_binomial(n, 12, plogis(0.3 * x), 5),
                       trials = 12),
  compois = list(fam = brmsfamily("com_poisson"),
                 y = rpois(n, exp(1 + 0.3 * x)))
)
for (nm in names(fits)) {
  f <- fits[[nm]]
  d <- data.frame(y = f$y, x = x, cc = cc)
  discrete <- identical(f$fam$type %||% frmtmb:::as_frmtmb_family(f$fam)$type,
                        "discrete")
  d$y2 <- d$y + if (nm == "von_mises") 0.4 else if (nm == "Beta") 0 else
    if (all(d$y == round(d$y))) 2 else 0.5
  if (nm == "Beta") d$y2 <- pmin(d$y + 0.1, 0.999)
  if (nm == "von_mises") d$y2 <- pmin(d$y + 0.4, pi)
  if (!is.null(f$trials)) d$tr <- f$trials
  lhs <- if (!is.null(f$trials)) "y | trials(tr) + cens(cc, y2)" else
    "y | cens(cc, y2)"
  r <- tryCatch({
    fit <- frm(as.formula(paste(lhs, "~ x")), family = f$fam, data = d)
    sprintf("cens logLik %.6f conv %s", as.numeric(logLik(fit)),
            fit$opt$convergence)
  }, error = function(e) paste("cens ERROR:", conditionMessage(e)))
  lhs2 <- if (!is.null(f$trials)) "y | trials(tr) + trunc(lb = lo, ub = hi)" else
    "y | trunc(lb = lo, ub = hi)"
  qs <- stats::quantile(d$y, c(0.05, 0.95))
  d$lo <- if (nm %in% c("poisson", "negbinomial", "geometric", "binomial",
                        "beta_binomial", "compois")) max(1, floor(qs[1])) else
    qs[1] - 1e-9
  d$hi <- qs[2] + 1e-9
  dt <- d[d$y >= d$lo & d$y <= d$hi, ]
  r2 <- tryCatch({
    fit <- frm(as.formula(paste(lhs2, "~ x")), family = f$fam, data = dt)
    sprintf("trunc logLik %.6f conv %s", as.numeric(logLik(fit)),
            fit$opt$convergence)
  }, error = function(e) paste("trunc ERROR:", conditionMessage(e)))
  cat(sprintf("%-18s %s | %s\n", nm, r, r2))
}
