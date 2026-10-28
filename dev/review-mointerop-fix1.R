## Reviewer: short check of fix round 1.
##   Rscript dev/review-mointerop-fix1.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
q <- function(e) suppressWarnings(suppressMessages(e))
se_of <- function(fit, nd, ...) predictions(fit, newdata = nd, ...)$std.error

## 1a. gap fit (seed 7, the review's data, g in 20 groups of 10)
set.seed(7)
dg <- data.frame(x = sample(c(1, 2, 5), 200, TRUE))
dg$y <- 1 + 0.8 * dg$x + rnorm(200)
dg$g <- factor(rep(1:20, each = 10))
fg <- q(frm(y ~ mo(x), data = dg))
mg <- lm(y ~ factor(x), data = dg)
for (v in list("HC0", "HC1", "stata", ~ g)) {
  lab <- if (inherits(v, "formula")) "~ g" else v
  show(paste("gap", lab, "frm x = 1..5 / lm x = 1, 2, 5"), rbind(
    frm = se_of(fg, data.frame(x = 1:5), vcov = v),
    lm = c(se_of(mg, data.frame(x = c(1, 2)), vcov = v), NA, NA,
           se_of(mg, data.frame(x = 5), vcov = v))))
}
show("gap full matrix", se_of(fg, data.frame(x = c(1, 2, 5)),
                              vcov = vcov(fg, full = TRUE)))

## 1b. the conditional bread against a direct computation: a generalized
## inverse of the outer Hessian. For quantities that do not move along
## the null direction, J ginv(H) S'S ginv(H) J' does not depend on the
## generalized inverse chosen.
frm_ns <- asNamespace("frmtmb")
show("direct conditional-bread sandwich vs get_vcov(HC0) (gap fit)", {
  om <- frm_ns$outer_par_map(fg)
  par <- fg$opt$par
  H <- numDeriv::jacobian(fg$obj$gr, par)
  H <- (H + t(H)) / 2
  S <- frm_ns$cluster_scores_at(fg, frm_ns$resolve_cluster(
    fg, seq_len(200)))
  Hi <- MASS::ginv(H, tol = 1e-7)
  Vg <- Hi %*% crossprod(S) %*% Hi
  mi <- fg$frame$linpreds[[1]]$mo[[1]]
  gaps <- frm_ns$interop_mo_gaps(fg)[[mi$zeta]]
  fmap <- function(p) {
    z <- p[om$comp == mi$zeta]
    c(p[om$comp %in% c("beta", "betad")],
      cumsum(frm_ns$mo_simplex(z, mi$chart))[gaps$obs])
  }
  J <- numDeriv::jacobian(fmap, par)
  Vd <- J %*% Vg %*% t(J)
  Vl <- get_vcov(fg, vcov = "HC0")
  list(eig_H = signif(eigen(H, only.values = TRUE)$values, 4),
       direct_se = sqrt(diag(Vd)), lane_se = sqrt(diag(Vl)),
       max_rel = max(abs(Vd - Vl)) / max(abs(Vl)))
})

## 1c. one-value fit (seed 5)
set.seed(5)
n <- 200
d1 <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                 g = factor(sample(1:8, n, TRUE)))
d1$y <- 1 + 0.8 * d1$x + 0.5 * d1$z + rnorm(8)[d1$g] * 0.5 + rnorm(n)
d1$x1 <- 4L
f1 <- q(frm(y ~ mo(x1) + z, data = d1))
m1 <- lm(y ~ z, data = d1)
nd1 <- data.frame(x1 = 4L, z = c(0, 1))
for (v in list("HC0", "HC1", ~ g)) {
  lab <- if (inherits(v, "formula")) "~ g" else v
  show(paste("one-value", lab), rbind(frm = se_of(f1, nd1, vcov = v),
                                      lm = se_of(m1, nd1["z"], vcov = v)))
}
show("one-value full matrix", se_of(f1, nd1, vcov = vcov(f1, full = TRUE)))

## 1d. boundary RE fit (seed 31 of -vcov.R)
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:30, n, TRUE)),
                w = runif(n, 0.5, 2), off = log(runif(n, 1, 3)))
d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x))
d$ybin <- rbinom(n, 1, pnorm(0.2 + 0.6 * d$x))
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + d$off))
d$y <- 1 + 0.5 * d$x + rnorm(n) * (1 + abs(d$x))
fr <- q(frm(y ~ x + (1 | g), data = d))
show("boundary RE ~ g vs glm(gaussian) vcovCL HC1 / CR1 ratio", {
  a <- sqrt(diag(get_vcov(fr, vcov = ~ g)))[1:2]
  m <- glm(y ~ x, data = d)
  b <- sqrt(diag(sandwich::vcovCL(m, cluster = ~ g)))
  list(frm = a, glm = b, ratio = a / b)
})
show("boundary RE HC0 (refused: RE crosses rows)", get_vcov(fr, vcov = "HC0"))
show("boundary RE HC1", get_vcov(fr, vcov = "HC1"))

## 1e. non-gap fits: sandwich vs rellib-r14's vcov_cluster (see -fix1b.R)
f0 <- frm(y ~ x, data = d)
fp <- frm(cnt ~ x + offset(off), data = d, family = poisson())
dm <- data.frame(k = sample(1:4, n, TRUE), z = rnorm(n),
                 g = factor(sample(1:30, n, TRUE)))
dm$y <- 1 + c(0, 0.3, 0.5, 0.9)[dm$k] + 0.4 * dm$z + rnorm(n)
fm <- frm(y ~ mo(k) + z, data = dm)
set.seed(33)
dn <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
dn$y <- 1 + 0.5 * dn$x + rnorm(30)[dn$g] + rnorm(n)
fn <- frm(y ~ x + (1 | g), data = dn)
out <- list(
  f0_hc0 = vcov_cluster(f0, seq_len(n), "CR0", full = TRUE),
  f0_g = vcov_cluster(f0, ~ g, "CR1", full = TRUE),
  fp_g = vcov_cluster(fp, ~ g, "CR1", full = TRUE),
  fm_g = vcov_cluster(fm, ~ g, "CR1", full = TRUE),
  fn_g = vcov_cluster(fn, ~ g, "CR1", full = TRUE),
  f0_me = get_vcov(f0, vcov = "HC0"), fp_me = get_vcov(fp, vcov = ~ g),
  fn_me = get_vcov(fn, vcov = ~ g))
saveRDS(out, "C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-out/fix1-nongap-lane.rds")

## 2. one-value: df, AIC, anova, predict(), (mo(x1) | g)
show("df / AIC frm vs lm", rbind(
  frm = c(df = attr(logLik(f1), "df"), AIC = AIC(f1), BIC = BIC(f1)),
  lm = c(df = attr(logLik(m1), "df"), AIC = AIC(m1), BIC = BIC(m1))))
show("anova", anova(frm(y ~ z, data = d1), f1))
show("predict() one-value", {
  set.seed(1)
  a <- predict(f1, newdata = nd1, ndraws = 4000)
  set.seed(1)
  b <- predict(frm(y ~ z, data = d1), newdata = nd1["z"], ndraws = 4000)
  rbind(one_value = a[, "Est.Error"], y_z = b[, "Est.Error"])
})
show("(mo(x1) | g)", frm(y ~ mo(x1) + z + (mo(x1) | g), data = d1))
show("(1 | g) with mo(x1) is fine", fixef(q(frm(y ~ mo(x1) + z + (1 | g),
                                                data = d1))))

## 3. HC1 k, REML refusal, lme4/glmmTMB behavior
show("HC1 k: gap fit, one-value fit, y ~ x, mo(k) + z", c(
  gap = frm_ns$interop_hc1_k(fg), one = frm_ns$interop_hc1_k(f1),
  lin = frm_ns$interop_hc1_k(f0), mo4 = frm_ns$interop_hc1_k(fm)))
show("HC1/HC0 on mo(k) + z vs lm(factor(k) + z)", rbind(
  frm = sqrt(diag(get_vcov(fm, vcov = "HC1")) /
               diag(get_vcov(fm, vcov = "HC0")))[1],
  lm = sqrt(diag(sandwich::vcovHC(lm(y ~ factor(k) + z, data = dm),
                                  "HC1")) /
              diag(sandwich::vcovHC(lm(y ~ factor(k) + z, data = dm),
                                    "HC0")))[1]))
show("sigma ~ x: HC1 k counts beta only", {
  fs <- frm(bf(y ~ x, sigma ~ x), data = d)
  c(k = frm_ns$interop_hc1_k(fs), n = n)
})
frl <- q(frm(y ~ x + (1 | g), data = dn, REML = TRUE))
show("REML full matrix", get_vcov(frl, vcov = vcov(frl, full = TRUE)))
show("REML 2 * get_vcov", sqrt(diag(get_vcov(frl, vcov = 2 * get_vcov(frl)))))
ml4 <- lme4::lmer(y ~ x + (1 | g), data = dn, REML = FALSE)
show("lme4 + marginaleffects HC1", predictions(
  ml4, newdata = data.frame(x = 0, g = factor(1, levels = 1:30)),
  vcov = "HC1")$std.error)
show("lme4 + marginaleffects ~ g", predictions(
  ml4, newdata = data.frame(x = 0, g = factor(1, levels = 1:30)),
  vcov = ~ g)$std.error)
mt <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = dn)
show("glmmTMB + marginaleffects HC1", predictions(
  mt, newdata = data.frame(x = 0, g = factor(1, levels = 1:30)),
  vcov = "HC1")$std.error)
show("glmmTMB + marginaleffects HC0", predictions(
  mt, newdata = data.frame(x = 0, g = factor(1, levels = 1:30)),
  vcov = "HC0")$std.error)
show("frm RE HC1", get_vcov(fn, vcov = "HC1"))

## 5. get_predict refusals and pass-through
fs <- frm(bf(y ~ x, sigma ~ x), data = d)
nds <- data.frame(x = c(-1, 0, 1))
show("dpar = sigma", rbind(
  me = q(predictions(fs, newdata = nds, dpar = "sigma"))$estimate,
  fitted = fitted(fs, newdata = nds, dpar = "sigma")[, 1]))
show("dpar = sigma SE", rbind(
  me = q(predictions(fs, newdata = nds, dpar = "sigma"))$std.error,
  fitted = fitted(fs, newdata = nds, dpar = "sigma")[, 2]))
show("avg_slopes dpar = sigma", q(avg_slopes(fs, variables = "x",
                                              dpar = "sigma")))
set.seed(3)
x <- rnorm(n)
dmv <- data.frame(x = x, y1 = 1 + x + rnorm(n), y2 = -1 + 0.5 * x + rnorm(n))
fmv <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)), data = dmv)
show("resp = y2", rbind(
  me = q(predictions(fmv, newdata = nds, resp = "y2"))$estimate,
  fitted = fitted(fmv, newdata = nds, resp = "y2")[, 1]))
show("no resp on mv fit", head(as.data.frame(q(predictions(
  fmv, newdata = nds)))[, c("rowid", "group", "estimate")]))
show("re.form refused", predictions(fn, newdata = data.frame(
  x = 0, g = factor(1, levels = 1:30)), re.form = NA))
show("ndraws refused", predictions(fn, newdata = data.frame(
  x = 0, g = factor(1, levels = 1:30)), ndraws = 10))
show("ordinary calls still run: avg_comparisons by, hypotheses, wts", {
  dn$w <- runif(n)
  list(avg_comparisons(fn, variables = "x", by = "g")[1:2, "estimate"],
       hypotheses(fn, "b2 = 0")$estimate,
       avg_predictions(fn, wts = "w")$estimate)
})

## 6. seed-23 fallback and 7. warning count
set.seed(23)
n <- 300
d3 <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                 g = factor(sample(1:10, n, TRUE)))
d3$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d3$g]) * d3$x + 0.5 * d3$z + rnorm(n)
f23 <- q(frm(y ~ mo(x) + z + (mo(x) | g), data = d3))
show("seed 23 get_vcov", get_vcov(f23))
show("seed 23 predictions observed rows", se_of(f23, data.frame(
  x = c(1, 2, 5), z = 0, g = factor(1, levels = 1:10))))
set.seed(51)
n <- 600
G <- 20
d5 <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                 g = factor(sample(seq_len(G), n, TRUE)))
d5$y <- 1 + rnorm(G, 0, 0.6)[d5$g] + 0.4 * (d5$x - 1) + 0.5 * d5$z +
  rnorm(n)
f51 <- q(frm(y ~ mo(x) + z + (1 | g), data = d5))
show("warning count, 4 rows, 20-level g", se_of(f51, data.frame(
  x = c(1, 2, 3, 5), z = 0, g = factor(1, levels = 1:G))))
