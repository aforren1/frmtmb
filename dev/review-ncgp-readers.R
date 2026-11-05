# Item 5: readers on both builds at ordinary noise (sigma 0.3). The answer should
# not move above the nugget's own effect. Writes one RDS per build; compared by
# review-ncgp-readers-cmp.R. Args: lane|base
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs(arm)
set.seed(8150)
n <- 90
d <- data.frame(x = runif(n, 0, 6), z = runif(n, 0, 6), g = factor(sample(letters[1:6], n, TRUE)),
                f = factor(sample(c("u", "v", "w"), n, TRUE)), w = runif(n, 0.5, 2), xr = rep(sort(runif(18, 0, 6)), each = 5))
d$y <- sin(d$x) + 0.3 * d$z + rnorm(6)[as.integer(d$g)] * 0.4 + rnorm(n, sd = 0.3)
d$y2 <- cos(d$x) + rnorm(n, sd = 0.3)
d$cnt <- rpois(n, exp(1 + 0.7 * sin(d$x))); d$bin <- rbinom(n, 1, plogis(sin(d$x)))
d$yr <- sin(d$xr) + rnorm(n, sd = 0.3)
d$tr <- 20; d$bn <- rbinom(n, 20, plogis(sin(d$x)))
nd <- data.frame(x = seq(0.25, 5.75, length.out = 9) + 0.013, z = 3, g = "a", f = "u", w = 1, xr = 3.21)
specs <- list(
  gp = list(f = y ~ gp(x)), gpz = list(f = y ~ z + gp(x)), re = list(f = y ~ z + gp(x) + (1 | g)),
  reml = list(f = y ~ z + gp(x), REML = TRUE), by = list(f = y ~ gp(x, by = f)), byw = list(f = y ~ gp(x, by = w)),
  grF = list(f = yr ~ gp(xr, gr = FALSE)), rep = list(f = yr ~ gp(xr)), iso2 = list(f = y ~ gp(x, z)), aniso = list(f = y ~ gp(x, z, iso = FALSE)),
  cmcF = list(f = y ~ gp(x, by = f, cmc = FALSE)), sigma = list(f = bf(y ~ z, sigma ~ gp(x))),
  nl = list(f = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)), mv = list(f = bf(y ~ gp(x)) + bf(y2 ~ gp(x)) + set_rescor(FALSE)),
  two = list(f = y ~ gp(x) + gp(z)), pois = list(f = cnt ~ gp(x), family = poisson()), bern = list(f = bin ~ gp(x), family = bernoulli()),
  binom = list(f = bn | trials(tr) ~ gp(x), family = binomial()), student = list(f = y ~ gp(x), family = student()), lognorm = list(f = exp(y / 3) ~ gp(x), family = lognormal())
)
tr <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(err) paste("ERR:", substr(conditionMessage(err), 1, 110)))
out <- list()
for (nm in names(specs)) {
  s <- specs[[nm]]; a <- s; a$f <- NULL
  fit <- tr(do.call(frm, c(list(s$f, data = d), a)))
  if (is.character(fit)) { out[[nm]] <- list(fit = fit); cat(nm, fit, "\n"); next }
  var <- if (grepl("xr", deparse(s$f)[1])) "xr" else "x"
  rd <- list(
    ll = as.numeric(logLik(fit)), conv = fit$opt$convergence,
    fixef = tr(fixef(fit, flatten = TRUE)), vc = tr(unlist(VarCorr(fit))), se = tr(sqrt(diag(as.matrix(vcov(fit))))),
    theta = tr(unname(fit$estimates$theta)),
    ci = tr(confint(fit)),
    fitted = tr(fitted(fit)), fittedn = tr(fitted(fit, newdata = nd)),
    linp = tr(frm_linpred(fit, newdata = nd, se.fit = TRUE)),
    ran = tr({ r <- ranef(fit); if (is.list(r)) unlist(r) else r }),
    coefs = tr(unlist(coef(fit))),
    resid = tr(residuals(fit)), simmean = tr({ set.seed(1); colMeans(as.matrix(simulate(fit, nsim = 200, seed = 1))) }),
    curve = tr({ cv <- frm_curve(fit, var = var, newdata = nd); as.data.frame(cv) }),
    deriv = tr({ cv <- frm_curve_deriv(fit, var = var, newdata = nd); as.data.frame(cv) }),
    deriv2 = tr({ cv <- frm_curve_deriv(fit, var = var, order = 2L, newdata = nd); as.data.frame(cv) }),
    ce = tr({ ce <- conditional_effects(fit, effects = var, plot = FALSE); as.data.frame(ce[[1]])[, c("estimate__", "se__", "lower__", "upper__")] }),
    emm = tr({ as.data.frame(emmeans::emmeans(fit, ~ 1))[, 1:3] }),
    emtr = tr({ as.data.frame(emmeans::emtrends(fit, ~ 1, var = var))[, 1:3] }),
    mfx = tr({ as.data.frame(marginaleffects::avg_slopes(fit, variables = var))[, c("estimate", "std.error")] }),
    mpred = tr({ as.data.frame(marginaleffects::predictions(fit, newdata = nd))[, c("estimate", "std.error")] }),
    anv = tr({ a0 <- anova(fit); as.numeric(unlist(a0)) }),
    secheck = tr(unlist(diagnose(fit, print = FALSE)$conv %||% 0))
  )
  out[[nm]] <- rd
  cat(nm, "done ll", rd$ll, "conv", rd$conv, "\n"); flush.console()
}
saveRDS(out, file.path(REV_OUT, sprintf("readers-%s.rds", arm)))
cat("DONE\n")
