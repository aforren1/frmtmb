# Attack 5: every consumer of RE structure and of the fit, on latent fits.
# Each call is reported as ok, a refusal, or an error, with what it printed.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
mk <- function() {
  set.seed(167)
  d <- expand.grid(t = 1:6, g = factor(1:25))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  d$z <- rnorm(n)
  d$f <- factor(sample(letters[1:3], n, TRUE))
  e <- as.vector(replicate(25, arima.sim(list(ar = 0.5), 6, sd = 0.5)))
  re <- rep(rnorm(25, 0, .4), each = 6)
  d$y <- rpois(n, exp(0.5 + 0.4 * d$x + e + re))
  d$bin <- rbinom(n, 1, plogis(0.3 * d$x + e))
  d$pos <- rgamma(n, 4, 4 / exp(0.5 + 0.4 * d$x + e))
  d
}
d <- mk()
fits <- list(
  arT = frm(y ~ x + z + ar(t, g, cov = TRUE), family = poisson(), data = d),
  reAR = frm(y ~ x + (1 | g) + ar(t, g, cov = TRUE), family = poisson(), data = d),
  slope = frm(y ~ x + z + (x | g) + ma(t, g, cov = TRUE), family = poisson(), data = d),
  arF = frm(y ~ x + ar(t, g), family = poisson(), data = d),
  gamF = frm(pos ~ x + ar(t, g), family = Gamma(link = "log"), data = d),
  cosy = frm(y ~ x + cosy(t, g), family = poisson(), data = d),
  unstr = frm(y ~ x + unstr(t, g), family = poisson(), data = d))

try1 <- function(label, expr) {
  r <- msg_of(eval.parent(substitute(expr)))
  v <- r$val
  s <- if (!is.null(r$err)) {
    paste0("ERR: ", substr(r$err, 1, 240))
  } else {
    txt <- tryCatch(paste(utils::capture.output(utils::str(v, max.level = 1, give.attr = FALSE, list.len = 4))[1:min(3, length(utils::capture.output(utils::str(v, max.level = 1, give.attr = FALSE, list.len = 4))))], collapse = " / "), error = function(e) "?")
    paste0("ok: ", substr(txt, 1, 200))
  }
  if (length(r$warn)) s <- paste0(s, "  [WARN: ", substr(paste(r$warn, collapse = " | "), 1, 200), "]")
  cat(sprintf("%-28s %s\n", label, s))
  invisible(v)
}

for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n=====", nm, " code", f$opt$convergence, " logLik", as.numeric(logLik(f)), "\n")
  try1("print", capture.output(print(f)))
  try1("summary", capture.output(print(summary(f))))
  try1("fixef", fixef(f))
  r <- try1("ranef", ranef(f))
  cat("   ranef names:", names(r), "\n")
  r <- try1("coef", coef(f))
  cat("   coef names:", names(r), "\n")
  try1("VarCorr", VarCorr(f))
  try1("ngrps", ngrps(f))
  try1("nobs", nobs(f))
  try1("logLik", logLik(f))
  try1("AIC", AIC(f))
  try1("BIC", BIC(f))
  try1("vcov", vcov(f))
  try1("vcov full", vcov(f, full = TRUE))
  try1("confint", confint(f))
  try1("confint profile", confint(f, parm = "x", method = "profile"))
  try1("confint_varcorr", confint_varcorr(f))
  try1("profile", profile(f, which = "x"))
  try1("fitted", fitted(f))
  try1("fitted se", fitted(f, se.fit = TRUE))
  try1("fitted re NA", fitted(f, re_formula = NA))
  try1("fitted re ~1|g", fitted(f, re_formula = ~ (1 | g)))
  try1("linpred", frm_linpred(f))
  try1("linpred se", frm_linpred(f, se.fit = TRUE))
  try1("predict", predict(f, ndraws = 50))
  try1("residuals", residuals(f))
  try1("residuals pearson", residuals(f, type = "pearson"))
  try1("residuals osa", residuals(f, type = "osa"))
  try1("dharma", dharma_residuals(f))
  try1("simulate", simulate(f, nsim = 2))
  try1("simulate re NA", simulate(f, nsim = 2, re_formula = NA))
  nd <- d[d$g %in% c("1", "2"), ]
  try1("fitted newdata", fitted(f, newdata = nd))
  try1("fitted newdata NA", fitted(f, newdata = nd, re_formula = NA))
  try1("linpred newdata", frm_linpred(f, newdata = nd))
  try1("predict newdata", predict(f, newdata = nd, ndraws = 50))
  try1("simulate newdata", simulate(f, nsim = 2, newdata = nd))
  nd2 <- nd
  nd2$g <- factor("new")
  nd2$t <- seq_len(nrow(nd2))
  try1("predict new group", predict(f, newdata = nd2, ndraws = 50, allow_new_levels = TRUE))
  try1("fitted new group", fitted(f, newdata = nd2, allow_new_levels = TRUE))
  try1("update", update(f, . ~ . + f))
  try1("update newdata", update(f, newdata = d[d$t < 6, ]))
  try1("refit", refit(f))
  try1("drop1", drop1(f))
  try1("add1", add1(f, ~ . + f))
  try1("anova", anova(f, f))
  try1("step", step(f, trace = 0))
  try1("conditional_effects", conditional_effects(f, "x"))
  try1("emmeans", emmeans::emmeans(f, ~ x, at = list(x = 0)))
  try1("marginaleffects", marginaleffects::avg_slopes(f, variables = "x"))
  try1("getME Z", frmtmb::getME(f, "Z"))
  try1("getME theta", frmtmb::getME(f, "theta"))
  try1("getME b", frmtmb::getME(f, "b"))
  try1("variables", variables(f))
  try1("hypothesis", hypothesis(f, "sderr > 0.1"))
  try1("hypothesis ar", hypothesis(f, "x > 0"))
  try1("default_prior", default_prior(f))
  try1("prior_summary", prior_summary(f))
  try1("autocor_matrix", autocor_matrix(f))
  try1("diagnose", capture.output(print(diagnose(f))))
  try1("influence", cooks.distance(f))
  try1("hatvalues", hatvalues(f))
  try1("weights", weights(f))
  try1("model.frame", model.frame(f))
  try1("model.matrix", model.matrix(f))
  try1("formula", formula(f))
  try1("terms", terms(f))
  tmp <- tempfile(fileext = ".rds")
  saveRDS(f, tmp)
  f2 <- readRDS(tmp)
  try1("saved logLik", logLik(f2))
  try1("saved fitted", fitted(f2))
  try1("saved predict nd", predict(f2, newdata = nd, ndraws = 10))
  try1("saved simulate", simulate(f2, nsim = 1))
  try1("saved update", update(f2, . ~ . + f))
  cat("  size MB", file.size(tmp) / 1e6, "\n")
}

cat("\n===== quadrature and importance\n")
try1("quadrature ar", frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d, quadrature = 5))
try1("quadrature ar + RE", frm(y ~ x + (1 | g) + ar(t, g, cov = TRUE), family = poisson(), data = d, quadrature = 5))
try1("importance ar", frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d, importance = 100))
try1("REML ar", frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d, REML = TRUE))
try1("subset", frm(y | subset(t < 6) ~ x + ar(t, g, cov = TRUE), family = poisson(), data = transform(d, keep = t < 6)))
dna <- d
dna$x[c(3, 40, 41)] <- NA
fna <- try1("NA x omit", frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = dna))
if (!is.null(fna)) {
  cat("  length fitted", length(fitted(fna)[, 1]), "nobs", nobs(fna), "\n")
  try1("NA predict", predict(fna, ndraws = 10))
  try1("NA fitted exclude", {
    fx <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = dna, na.action = na.exclude)
    dim(fitted(fx))
  })
}
try1("factor time", frm(y ~ x + ar(ft, g, cov = TRUE), family = poisson(), data = transform(d, ft = factor(t))))
try1("two autocor", frm(y ~ x + ar(t, g, cov = TRUE) + cosy(t, g), family = poisson(), data = d))
try1("dpar term", frm(bf(y ~ x, shape ~ ar(t, g, cov = TRUE)), family = negbinomial(), data = d))
try1("spec car etc", frm(y ~ x + s(z) + ar(t, g, cov = TRUE), family = poisson(), data = d))
try1("gr() ar1", frm(y ~ x + ar1(factor(t) + 0 | g) + ar(t, g, cov = TRUE), family = poisson(), data = d))
try1("mv", frm(bf(y ~ x + ar(t, g, cov = TRUE)) + bf(bin ~ x, family = bernoulli()), family = poisson(), data = d))
cat("done\n")
