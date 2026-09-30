# Reviewer re-check (punch round 1), lane ordinal: B3's recycle0 = TRUE
# edits sit on shared naming paths. Every name a fit exposes, on a broad
# set of NON-ordinal models (and a few ordinal flexible ones), for an
# identical() comparison of base and lane. Includes zero-length blocks: no
# random effects, no dpar formula, a smooth whose null space is empty,
# intercept-only, no fixed effects at all (y ~ 0 + (1 | g)).
# Usage: Rscript ... <base|lane>. Data seed 20261014.
# Output: dev/ordinal-rev2-names-<arm>.rds
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261014)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), t = seq_len(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                h = factor(sample(c("p", "q", "r"), n, TRUE)))
re <- rnorm(8, 0, 0.5)[as.integer(d$g)]
d$yg <- 1 + 0.5 * d$x + sin(2 * d$z) + re + rnorm(n)
d$yg2 <- 0.3 * d$x + rnorm(n)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x + re))
d$yb <- rbinom(n, 1, stats::plogis(d$x + re))
d$ycat <- factor(sample(c("A", "B", "C"), n, TRUE))
d$yz <- ifelse(runif(n) < 0.3, 0L, rpois(n, exp(0.5 + 0.3 * d$x)))
d$ymix <- ifelse(runif(n) < 0.4, rnorm(n, -2), rnorm(n, 2)) + 0.5 * d$x
d$y <- 1L + findInterval(stats::rlogis(n) + d$x, c(-1, 0, 1))
models <- list(
  lm_like = quote(frm(yg ~ x, data = d)),
  intercept_only = quote(frm(yg ~ 1, data = d)),
  no_fixed = quote(frm(yg ~ 0 + (1 | g), data = d)),
  re_int = quote(frm(yg ~ x + (1 | g), data = d)),
  re_slope_cor = quote(frm(yg ~ x + (1 + x | g), data = d)),
  re_two = quote(frm(yc ~ x + (1 | g) + (1 | h), family = poisson(), data = d)),
  smooth_tp = quote(frm(yg ~ s(z), data = d)),
  smooth_null_empty = quote(frm(yg ~ x + s(z, bs = "tp", m = 1), data = d)),
  smooth_re = quote(frm(yg ~ s(g, bs = "re"), data = d)),
  sigma_formula = quote(frm(bf(yg ~ x, sigma ~ z), data = d)),
  sigma_re = quote(frm(bf(yg ~ x, sigma ~ (1 | g)), data = d)),
  fixed_nu = quote(frm(bf(yg ~ x, nu = 4), family = student(), data = d)),
  bernoulli_re = quote(frm(yb ~ x + (1 | g), family = bernoulli(), data = d)),
  categorical = quote(frm(ycat ~ x, family = categorical(), data = d)),
  zip_zi = quote(frm(bf(yz ~ x, zi ~ x), family = zero_inflated_poisson(),
                     data = d)),
  mixture = quote(frm(ymix ~ x, family = mixture(gaussian(), gaussian()),
                      data = d)),
  mv_gauss = quote(frm(bf(yg ~ x) + bf(yg2 ~ x), data = d)),
  mv_gauss_rescor = quote(frm(bf(yg ~ x) + bf(yg2 ~ x) + set_rescor(TRUE),
                              data = d)),
  ar1 = quote(frm(yg ~ x + ar(time = t), data = d)),
  negbin = quote(frm(yc ~ x, family = negbinomial(), data = d)),
  ord_flex = quote(frm(y ~ x, family = cumulative(), data = d)),
  ord_flex_re = quote(frm(y ~ x + (1 | g), family = acat(), data = d)),
  ord_mv_gauss = quote(frm(bf(y ~ x) + bf(yg ~ x), data = d,
                           family = list(sratio(), gaussian())))
)
grab <- function(expr) tryCatch(suppressWarnings(expr), error = function(e)
  paste("ERROR:", conditionMessage(e)))
out <- list()
for (nm in names(models)) {
  cat("model", nm, "\n")
  fit <- grab(eval(models[[nm]]))
  if (is.character(fit)) { out[[nm]] <- list(fit = fit); next }
  tp <- grab(par_template(fit))
  out[[nm]] <- list(
    variables = grab(variables(fit)),
    coef_names = grab(names(coef(fit))),
    coef_inner = grab(lapply(coef(fit), names)),
    confint_rows = grab(rownames(confint(fit))),
    brms_labels = grab(brms_par_labels(fit)),
    flatten = grab(names(fixef(fit, flatten = TRUE))),
    fixef_rows = grab(rownames(fixef(fit))),
    vcov_rows = grab(rownames(vcov(fit))),
    template = grab(lapply(unclass(tp), names)),
    template_lengths = grab(lengths(unclass(tp))),
    default_prior = grab(as.data.frame(default_prior(fit))),
    hyp = grab(rownames(hypothesis(fit, "Intercept = 0", class = "b")$hypothesis)),
    summary = grab(utils::capture.output(print(summary(fit)))))
  if (nm %in% c("re_int", "sigma_formula", "categorical", "mv_gauss",
                "smooth_null_empty", "no_fixed")) {
    ds <- grab(suppressMessages(frm_sample(fit, chains = 1, iter = 100,
                                           refresh = 0, seed = 5)))
    out[[nm]]$draws_vars <- grab(variables(ds))
    out[[nm]]$draws_fixef <- grab(rownames(fixef(ds)))
  }
}
saveRDS(out, file.path(wt, paste0("dev/ordinal-rev2-names-", arm, ".rds")))
cat("done\n")
