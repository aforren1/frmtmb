# Every reader of an exact gp() fit runs on the lane: a table of what each returns
# against a reference where one exists. Seed 20261106, 60 rows unless stated.
# Output NCGPR lines: reader, spelling, ok/refused/ERROR, a number.
#   Rscript dev/ncgp-readers.R
LB <- commandArgs(TRUE)
.libPaths(c(if (!length(LB)) "C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(emmeans); library(marginaleffects)})
set.seed(20261106)
n <- 60L
d <- data.frame(x = round(runif(n, 0, 6), 3), z = round(runif(n, 0, 4), 2),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                w = runif(n, 0.5, 2), g = factor(rep(1:6, length.out = n)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(2 * d$x) + 0.2 * rnorm(n) + ifelse(d$fac == "B", 0.5, 0)
d$y2 <- cos(d$z) + 0.2 * rnorm(n)
d$cnt <- rpois(n, exp(1 + sin(2 * d$x)))
d$bin <- rbinom(n, 1, plogis(sin(2 * d$x)))
d$ysg <- rnorm(n, 0, exp(-1 + sin(d$x)))
d$ynl <- (2 + sin(d$z)) * exp(0.3 * d$w) + 0.2 * rnorm(n)
nd <- data.frame(x = c(0.111, 1.777, 3.333), z = c(1, 2, 3),
                 f = factor(c("a", "b", "c"), levels = levels(d$f)),
                 w = c(1, 1, 1), g = factor(c(1, 2, 3), levels = levels(d$g)),
                 fac = factor("A", levels = levels(d$fac)))
ok <- function(tag, sp, expr) {
  t0 <- Sys.time()
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, substr(conditionMessage(c), 1, 60)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "ncgp_err"))
  el <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(r, "ncgp_err")) {
    cat(sprintf("NCGPR %-20s %-26s ERROR/refused: %s\n", tag, sp, substr(unclass(r), 1, 120)))
  } else {
    nm <- if (is.numeric(r) && length(r)) sprintf("%.5g (n=%d, NaN %d)", mean(abs(r), na.rm = TRUE), length(r), sum(is.nan(r))) else
      paste(class(r)[1], paste(dim(r), collapse = "x"))
    cat(sprintf("NCGPR %-20s %-26s ok %s %.1fs%s\n", tag, sp, nm, el,
                if (length(w)) paste0(" warn: ", paste(unique(w), collapse = " | ")) else ""))
  }
  invisible(r)
}
F <- list(
  "gp(x)" = bf(y ~ fac + gp(x)),
  "gp(x, z)" = bf(y ~ gp(x, z)),
  "iso = FALSE" = bf(y ~ gp(x, z, iso = FALSE)),
  "by = f" = bf(y ~ gp(x, by = f)),
  "by = f, cmc = FALSE" = bf(y ~ gp(x, by = f, cmc = FALSE)),
  "by = w" = bf(y ~ gp(x, by = w)),
  "gr = FALSE" = bf(y ~ gp(x, gr = FALSE)),
  "scale = FALSE" = bf(y ~ gp(x, scale = FALSE)),
  "two gp" = bf(y ~ gp(x) + gp(z)),
  "gp + (1 | g)" = bf(y ~ gp(x) + (1 | g)),
  "sigma ~ gp(x)" = bf(ysg ~ 1, sigma ~ gp(x)),
  "nl b1 ~ gp(z)" = bf(ynl ~ b1 * exp(b2 * w), b1 ~ gp(z), b2 ~ 1, nl = TRUE),
  "mvbf" = mvbf(bf(y ~ gp(x)), bf(y2 ~ gp(z)), rescor = FALSE))
extra <- list(
  poisson = list(bf(cnt ~ gp(x)), poisson()),
  bernoulli = list(bf(bin ~ gp(x)), bernoulli()),
  student = list(bf(y ~ gp(x)), student()),
  reml = list(bf(y ~ fac + gp(x)), gaussian(), REML = TRUE),
  map = list(bf(y ~ gp(x)), gaussian(),
             prior = set_prior("normal(0, 1)", class = "sdgp")))
fits <- list()
for (nm in names(F)) {
  fits[[nm]] <- ok("frm", nm, frm(F[[nm]], data = d, family = gaussian()))
}
for (nm in names(extra)) {
  e <- extra[[nm]]
  args <- c(list(e[[1]], data = d, family = e[[2]]), e[-(1:2)])
  fits[[nm]] <- ok("frm", nm, do.call(frm, args))
}
fits <- Filter(function(f) inherits(f, "frmtmb_fit"), fits)
for (nm in names(fits)) {
  fit <- fits[[nm]]
  mv <- nm == "mvbf"
  ok("fitted", nm, as.vector(fitted(fit)))
  ok("predict", nm, as.vector(predict(fit)))
  ok("frm_linpred se", nm, if (mv) frm_linpred(fit, resp = "y", se.fit = TRUE)$se.fit else frm_linpred(fit, se.fit = TRUE)$se.fit)
  ok("frm_linpred new", nm, if (mv) frm_linpred(fit, newdata = nd, resp = "y", se.fit = TRUE)$se.fit else frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit)
  ok("fitted new", nm, as.vector(fitted(fit, newdata = nd)))
  ok("predict new", nm, as.vector(predict(fit, newdata = nd)))
  ok("re_formula NA", nm, as.vector(fitted(fit, re_formula = NA)))
  ok("residuals", nm, as.vector(residuals(fit)))
  ok("simulate", nm, as.vector(as.matrix(simulate(fit, nsim = 2, seed = 1))))
  ok("simulate new", nm, as.vector(as.matrix(simulate(fit, nsim = 2, seed = 1, newdata = nd))))
  ok("ranef", nm, as.vector(ranef(fit)[[1L]]))
  ok("ranef condSD", nm, as.vector(attr(ranef(fit, condVar = TRUE)[[1L]], "condSD")))
  ok("coef", nm, as.vector(coef(fit)[[1L]]))
  ok("VarCorr", nm, unlist(VarCorr(fit)))
  ok("vcov", nm, as.vector(vcov(fit)))
  ok("vcov full", nm, as.vector(vcov(fit, full = TRUE)))
  ok("summary", nm, summary(fit))
  ok("confint wald", nm, as.vector(as.matrix(confint(fit))))
  ok("confint_varcorr", nm, as.vector(as.matrix(confint_varcorr(fit)[, c("lwr", "upr")])))
  ok("logLik", nm, as.numeric(logLik(fit)))
  ok("AIC", nm, AIC(fit))
  ok("conditional_effects", nm, conditional_effects(fit, effects = "x"))
  ok("frm_curve", nm, frm_curve(fit, newdata = nd, simultaneous = FALSE)$.se)
  ok("frm_curve_deriv", nm, frm_curve_deriv(fit, var = "x", newdata = nd, simultaneous = FALSE)$.se)
  ok("frm_curve diff", nm, frm_curve(fit, newdata = nd, contrast = nd[3:1, ], simultaneous = FALSE)$.se)
  ok("update", nm, -update(fit, data = d)$opt$objective)
  ok("sigma", nm, as.numeric(sigma(fit)))
  ok("bayes_R2", nm, as.vector(bayes_R2(fit)))
  ok("emmeans", nm, emmeans::emmeans(fit, "fac", data = d))
  ok("avg_predictions", nm, as.data.frame(marginaleffects::avg_predictions(fit, newdata = nd))$std.error)
  ok("predictions", nm, as.data.frame(marginaleffects::predictions(fit, newdata = nd))$std.error)
  ok("diagnose", nm, diagnose(fit))
}
fit <- fits[["gp(x)"]]
ok("anova", "gp(x) vs + z", {
  f2 <- update(fit, . ~ . + z)
  anova(fit, f2)
})
ok("confint profile", "gp(x)", as.vector(as.matrix(confint(fit, method = "profile", parm = "sigma"))))
ok("refit own data", "gp(x)", -update(fit)$opt$objective)
ok("drop1", "gp(x)", drop1(fit))
ok("bootstrap", "gp(x)", frm_curve(fit, newdata = nd, band = "boot", nboot = 4, seed = 1))
ok("importance", "gp(x)", frm(bf(y ~ gp(x)), data = d, importance = 10))
ok("quadrature", "gp(x)", frm(bf(y ~ gp(x)), data = d, quadrature = TRUE))
ok("vcov_cluster", "gp(x)", as.vector(vcov_cluster(fit, cluster = d$g)))
ok("getME b", "gp(x)", as.vector(getME(fit, "b")))
ok("start list", "gp(x)", -frm(bf(y ~ gp(x)), data = d, start = list(theta = c(0, 0)))$opt$objective)
