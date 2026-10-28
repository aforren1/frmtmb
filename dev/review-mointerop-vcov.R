## Reviewer: marginaleffects vcov = on cases the lane did not test.
##   Rscript dev/review-mointerop-vcov.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
cmp <- function(fit, m, v) {
  a <- sqrt(diag(unname(get_vcov(fit, vcov = v))))
  b <- sqrt(diag(unname(get_vcov(m, vcov = v))))
  k <- min(length(a), length(b))
  rbind(frm = a[seq_len(k)], ref = b[seq_len(k)],
        ratio = a[seq_len(k)] / b[seq_len(k)])
}
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:30, n, TRUE)),
                w = runif(n, 0.5, 2), off = log(runif(n, 1, 3)))
d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x))
d$ybin <- rbinom(n, 1, pnorm(0.2 + 0.6 * d$x))
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + d$off))
d$y <- 1 + 0.5 * d$x + rnorm(n) * (1 + abs(d$x))

## Gamma(log): non-canonical link
fg <- frm(ypos ~ x, data = d, family = Gamma(link = "log"))
mg <- glm(ypos ~ x, data = d, family = Gamma(link = "log"))
show("Gamma log HC0", cmp(fg, mg, "HC0"))
show("Gamma log ~ g", cmp(fg, mg, ~ g))
## binomial probit
fb <- frm(ybin ~ x, data = d, family = bernoulli(link = "probit"))
mb <- glm(ybin ~ x, data = d, family = binomial(link = "probit"))
show("binomial probit HC0", cmp(fb, mb, "HC0"))
## Poisson with an offset
fo <- frm(cnt ~ x + offset(off), data = d, family = poisson())
mo_ <- glm(cnt ~ x + offset(off), data = d, family = poisson())
show("poisson offset HC0", cmp(fo, mo_, "HC0"))
## gaussian with weights (frmtmb weights() as brms's)
fw <- frm(y | weights(w) ~ x, data = d)
mw <- lm(y ~ x, data = d, weights = w)
show("gaussian weights coef", rbind(frm = fixef(fw)[, 1], lm = coef(mw)))
show("gaussian weights HC0", cmp(fw, mw, "HC0"))
show("gaussian weights model-based", rbind(
  frm = sqrt(diag(get_vcov(fw)))[1:2],
  lm_ml = sqrt(diag(vcov(mw))) * sqrt((n - 2) / n)))
## gaussian sigma ~ x (distributional)
fs <- frm(bf(y ~ x, sigma ~ x), data = d)
show("distributional HC0 get_vcov", sqrt(diag(get_vcov(fs, vcov = "HC0"))))
show("distributional vcov_cluster CR0 each row", sqrt(diag(vcov_cluster(
  fs, seq_len(n), "CR0"))))
## nonlinear
d$ynl <- 2 * exp(-0.5 * abs(d$x)) + rnorm(n, 0, 0.2)
fn <- show("nonlinear fit", frm(bf(ynl ~ a * exp(-b * abs(x)), a ~ 1, b ~ 1,
                                  nl = TRUE), data = d,
                               prior = NULL))
mn <- nls(ynl ~ a * exp(-b * abs(x)), data = d, start = list(a = 2, b = 0.5))
show("nonlinear HC0", {
  a <- sqrt(diag(get_vcov(fn, vcov = "HC0")))
  b <- sqrt(diag(sandwich::sandwich(mn)))
  list(frm = a, nls = b)
})
## mo() with 3 categories, HC0 vs lm on factor (same model)
set.seed(32)
dm <- data.frame(k = sample(1:3, n, TRUE))
dm$y <- 1 + c(0, 0.5, 0.7)[dm$k] + rnorm(n) * dm$k
fm <- frm(y ~ mo(k), data = dm)
mm <- lm(y ~ factor(k), data = dm)
show("mo(3) HC0 predictions vs lm factor HC0", {
  nd <- data.frame(k = 1:3)
  rbind(frm = predictions(fm, newdata = nd, vcov = "HC0")$std.error,
        lm = predictions(mm, newdata = nd, vcov = "HC0")$std.error)
})
## RE fit: whole outer matrix equals the default
fr <- frm(y ~ x + (1 | g), data = d)
show("RE fit: vcov(full=TRUE) vs default get_vcov", {
  A <- get_vcov(fr)
  B <- get_vcov(fr, vcov = vcov(fr, full = TRUE))
  list(names = rownames(A), max_rel = max(abs(A - B) / abs(A)))
})
show("RE fit: predictions default vs vcov(full=TRUE)", {
  nd <- data.frame(x = c(-1, 1), g = factor(1, levels = levels(d$g)))
  rbind(default = predictions(fr, newdata = nd)$std.error,
        full = predictions(fr, newdata = nd,
                           vcov = vcov(fr, full = TRUE))$std.error,
        double = predictions(fr, newdata = nd,
                             vcov = 2 * get_vcov(fr))$std.error)
})
show("RE fit ~ g (nested clusters)", sqrt(diag(get_vcov(fr, vcov = ~ g))))
show("RE fit ~ g direct vcov_cluster", sqrt(diag(vcov_cluster(fr, ~ g, "CR1"))))
## REML fit
frl <- frm(y ~ x + (1 | g), data = d, REML = TRUE)
show("REML get_vcov default names", rownames(get_vcov(frl)))
show("REML vcov = ~ g", get_vcov(frl, vcov = ~ g))
show("REML vcov(full = TRUE)", get_vcov(frl, vcov = vcov(frl, full = TRUE)))
## A fit with priors: vcov_cluster refuses
fp <- frm(y ~ x, data = d, prior = set_prior("normal(0, 1)", class = "b"))
show("prior fit HC0", get_vcov(fp, vcov = "HC0"))
## HC1 is HC0 * n / (n - k): no hat matrix needed
f0 <- frm(y ~ x, data = d)
m0 <- lm(y ~ x, data = d)
show("lm HC1/HC0 ratio", sqrt(diag(sandwich::vcovHC(m0, type = "HC1")) /
                               diag(sandwich::vcovHC(m0, type = "HC0"))))
show("frm HC1 refusal text", get_vcov(f0, vcov = "HC1"))
show("frm 'stata'", get_vcov(f0, vcov = "stata"))
show("frm 'HC'", get_vcov(f0, vcov = "HC"))
show("frm 'OPG'", get_vcov(f0, vcov = "OPG"))
show("frm lowercase 'hc0'", sqrt(diag(get_vcov(f0, vcov = "hc0"))))
show("lm lowercase 'hc0' via marginaleffects", tryCatch(
  predictions(m0, newdata = data.frame(x = 0), vcov = "hc0")$std.error,
  error = function(e) conditionMessage(e)))
## a gaussian fit: full unnamed matrix when outer and coef sizes agree
show("sizes", c(coef = length(get_coef(f0)),
                outer = nrow(vcov(f0, full = TRUE))))
show("unnamed full when sizes agree", {
  V <- unname(vcov(f0, full = TRUE))
  sqrt(diag(get_vcov(f0, vcov = V)))
})
show("named full", sqrt(diag(get_vcov(f0, vcov = vcov(f0, full = TRUE)))))
## comparisons with ~ g on a glm vs glm
show("avg_comparisons poisson ~ g vs glm", {
  rbind(frm = avg_slopes(fo, variables = "x", vcov = ~ g)$std.error,
        glm = avg_slopes(mo_, variables = "x", vcov = ~ g)$std.error)
})
## sigma reported as get_coef in the vcov?
show("gaussian get_coef names", names(get_coef(f0)))
