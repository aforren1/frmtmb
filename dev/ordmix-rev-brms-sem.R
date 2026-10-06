# Reviewer of lane ordmix: brms 2.23.0's own behavior on the lane's
# claims, by stancode(), standata(), default_prior() and brms internals.
# No compile. Output: dev/ordmix-rev-log-brms-sem.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sample(1:4, n, TRUE)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
d$w <- runif(n, 0.5, 2)
d$y2 <- sample(1:3, n, TRUE)
hdr <- function(s) cat("\n=====", s, "=====\n")
grab <- function(code, pat) {
  l <- strsplit(as.character(code), "\n")[[1]]
  cat(paste0("  | ", l[grepl(pat, l)]), sep = "\n")
}
try_sc <- function(f, fam, pat, ...) {
  r <- tryCatch(withCallingHandlers(
    stancode(f, data = d, family = fam, ...),
    message = function(m) {
      cat("  MESSAGE:", conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      cat("  WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }), error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    })
  if (!is.null(r)) grab(r, pat)
  invisible(r)
}
dp <- function(f, fam) {
  p <- tryCatch(suppressMessages(default_prior(f, data = d, family = fam)),
                error = function(e) {
                  cat("  ERROR:", conditionMessage(e), "\n")
                  NULL
                })
  if (!is.null(p)) {
    print(as.data.frame(p)[, c("prior", "class", "coef", "group", "dpar",
                               "lb")], row.names = FALSE)
  }
}

hdr("1a default order, two cumulative")
try_sc(bf(y ~ x), mixture(cumulative(), cumulative()),
       "Intercept|Xc|means_X|theta")
dp(bf(y ~ x), mixture(cumulative(), cumulative()))

hdr("1b order = mu, cumulative + sratio")
try_sc(bf(y ~ x), mixture(cumulative(), sratio(), order = "mu"),
       "Intercept|Xc|means_X")
dp(bf(y ~ x), mixture(cumulative(), sratio(), order = "mu"))

hdr("1c order = mu, sratio + acat (vector?)")
try_sc(bf(y ~ x), mixture(sratio(), acat(), order = "mu"), "Intercept")

hdr("1d order = mu, sum_to_zero component")
try_sc(bf(y ~ x), mixture(cumulative(), cratio(threshold = "sum_to_zero"),
                          order = "mu"), "Intercept|Xc")

hdr("1e order = mu with equidistant")
try_sc(bf(y ~ x), mixture(cumulative(threshold = "equidistant"),
                          cumulative(), order = "mu"), "delta")

hdr("1f equidistant under none (defect 4)")
try_sc(bf(y ~ x), mixture(cumulative(threshold = "equidistant"),
                          sratio(threshold = "equidistant")), "delta")
dp(bf(y ~ x), mixture(cumulative(threshold = "equidistant"),
                      sratio(threshold = "equidistant")))

hdr("1g disc ~ cs(z), sratio")
try_sc(bf(y ~ x, disc ~ cs(z)), sratio(), "cs")
hdr("1g2 disc1 ~ cs(z), mixture")
try_sc(bf(y ~ x, disc1 ~ cs(z)), mixture(sratio(), sratio()), "cs")
hdr("1g3 mu2 ~ cs(z) mixture sratio")
try_sc(bf(y ~ x, mu2 ~ x + cs(z)), mixture(cumulative(), sratio()),
       "cs")
hdr("1g4 sigma ~ cs(z), gaussian")
try_sc(bf(x ~ z, sigma ~ cs(y)), gaussian(), "cs")
hdr("1g5 cs() in a gaussian mu")
try_sc(bf(x ~ cs(y)), gaussian(), "cs")
hdr("1g6 cs() in an nlpar of sratio")
try_sc(bf(y ~ a, a ~ cs(z), nl = TRUE), sratio(), "cs")
hdr("1g7 cs() in hu of hurdle_cumulative")
try_sc(bf(yh ~ x, hu ~ cs(z)), hurdle_cumulative(), "cs")

hdr("1h hurdle thres(gr = )")
try_sc(bf(yh | thres(gr = g) ~ x), hurdle_cumulative(),
       "merged|nthres|Intercept")
dp(bf(yh | thres(gr = g) ~ x), hurdle_cumulative())
hdr("1i hurdle cs() with thres(gr = )")
try_sc(bf(yh | thres(gr = g) ~ cs(x)), hurdle_cumulative(), "cs")
hdr("1i2 hurdle cs() probit")
try_sc(bf(yh ~ cs(x)), hurdle_cumulative("probit"), "lpmf\\(Y|mucs")

hdr("1j mixture cumulative + hurdle_cumulative, y with zeros")
r <- try_sc(bf(yh ~ x), mixture(cumulative(), hurdle_cumulative()),
            "lpmf|nthres|Intercept")
sd1 <- tryCatch(standata(bf(yh ~ x), data = d,
                         family = mixture(cumulative(),
                                          hurdle_cumulative())),
                error = function(e) {
                  cat("  standata ERROR:", conditionMessage(e), "\n")
                  NULL
                })
if (!is.null(sd1)) cat("  standata nthres", sd1$nthres, "range Y",
                       range(sd1$Y), "\n")
hdr("1j2 same, y without zeros")
sd2 <- standata(bf(y ~ x), data = d,
                family = mixture(cumulative(), hurdle_cumulative()))
cat("  standata nthres", sd2$nthres, "range Y", range(sd2$Y), "\n")
hdr("1j3 hurdle_cumulative alone, y with zeros")
sd3 <- standata(bf(yh ~ x), data = d, family = hurdle_cumulative())
cat("  standata nthres", sd3$nthres, "range Y", range(sd3$Y), "\n")
hdr("1j4 cumulative alone, y with zeros")
tryCatch(standata(bf(yh ~ x), data = d, family = cumulative()),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

hdr("1k thres(gr = ) inside a mixture, unequal counts")
d$yu <- ifelse(d$g == "b", pmin(d$y, 3L), d$y)
try_sc(bf(yu | thres(gr = g) ~ x), mixture(cumulative(), cumulative()),
       "nthres|Intercept_mu1")
sdk <- standata(bf(yu | thres(gr = g) ~ x), data = d,
                family = mixture(cumulative(), cumulative()))
cat("  nthres", sdk$nthres, "\n")

hdr("1l weights, cens, trunc, mi on an ordinal mixture")
try_sc(bf(y | weights(w) ~ x), mixture(cumulative(), sratio()), "weights")
d$ce <- sample(c(0, 1), n, TRUE)
try_sc(bf(y | cens(ce) ~ x), mixture(cumulative(), sratio()), "cens")
try_sc(bf(y | cens(ce) ~ x), cumulative(), "cens")
try_sc(bf(y | trunc(ub = 4) ~ x), mixture(cumulative(), sratio()), "trunc")
try_sc(bf(y | trunc(ub = 4) ~ x), cumulative(), "trunc")

hdr("1m multivariate model with one ordinal-mixture response")
try_sc(mvbf(bf(y ~ x, family = mixture(cumulative(), cumulative())),
            bf(x ~ z, family = gaussian())), NULL, "lpmf|Intercept_y")

hdr("1n disc<k> default prior")
dp(bf(y ~ x, disc1 ~ z), mixture(cumulative(), cumulative()))
dp(bf(y ~ x, disc ~ z), cumulative())

hdr("1o is_ordinal / polytomous on a mixture")
mf <- mixture(cumulative(), cumulative())
cat("  is_ordinal", brms:::is_ordinal(mf), " is_polytomous",
    brms:::is_polytomous(mf), "\n")
cat("  incl_thres refusal text present:",
    any(grepl("not supported for mixture",
              deparse(brms:::prepare_predictions_thres))), "\n")

hdr("1p theta with predictors, three components")
try_sc(bf(y ~ x, theta1 ~ z), mixture(cumulative(), cumulative(),
                                      sratio()), "theta")
hdr("1q links line")
cat("  ", brms:::summarise_links(brmsterms(bf(y ~ x, disc1 ~ z)),
                                  wsp = 0), "\n")
