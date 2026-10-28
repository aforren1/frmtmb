## frmtmb on mo() with two values (D = 1) and one value (D = 0), on the
## data of dev/mointerop-brms-mo12.R (seed 5, 200 rows).
##   Rscript dev/mointerop-mo12.R <lane|base> > dev/mointerop-out/mo12-<arm>.txt
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
suppressMessages({
  library(frmtmb)
  library(marginaleffects)
  library(emmeans)
})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
show <- function(tag, expr) {
  cat("\n==", tag, "==\n")
  w <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(r)
  for (m in w) cat("WARNING:", substr(m, 1, 400), "\n")
  invisible(r)
}
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
d$xo <- factor(d$x, levels = 1:2, ordered = TRUE)
d$x1 <- 4
d$x01 <- d$x - 1L

f2 <- show("fit y ~ mo(x) + z", frm(y ~ mo(x) + z, data = d))
if (!inherits(f2, "error")) {
  ml <- lm(y ~ x01 + z, data = d)
  s <- summary(ml)
  sig_ml <- sqrt(sum(resid(ml)^2) / n)
  se_ml <- sqrt(diag(vcov(ml))) * sig_ml / s$sigma
  cat("\nexact (lm on the code, ML sigma):\n")
  print(cbind(est = coef(ml), se = se_ml), digits = 10)
  cat("logLik lm:", format(logLik(ml), digits = 12), "\n")
  show("fixef", fixef(f2))
  cat("logLik frm:", format(as.numeric(logLik(f2)), digits = 12), "df",
      attr(logLik(f2), "df"), "\n")
  show("summary", summary(f2))
  show("variables", variables(f2))
  show("default_prior", default_prior(y ~ mo(x) + z, data = d))
  show("prior_summary", prior_summary(f2))
  show("fitted newdata x = 1, 2, z = 0",
       fitted(f2, newdata = data.frame(x = c(1, 2), z = 0)))
  show("exact fitted SE at x = 1, 2, z = 0", predict(
    ml, newdata = data.frame(x01 = c(0, 1), z = 0), se.fit = TRUE)$se.fit *
      sig_ml / s$sigma)
  show("fitted x = 3", fitted(f2, newdata = data.frame(x = 3, z = 0)))
  show("fitted x = 0", fitted(f2, newdata = data.frame(x = 0, z = 0)))
  show("fitted x = 1.5", fitted(f2, newdata = data.frame(x = 1.5, z = 0)))
  ce <- show("conditional_effects x",
             conditional_effects(f2, effects = "x")$x[, c("x", "z",
                                                          "estimate__",
                                                          "se__")])
  show("conditional_effects z: x held at",
       unique(conditional_effects(f2, effects = "z")$z$x))
  show("conditional_effects z:x moderator values",
       unique(conditional_effects(f2, effects = "z:x")[[1]]$x))
  show("predictions x = 1, 2", as.data.frame(predictions(
    f2, newdata = data.frame(x = c(1, 2), z = 0)))[, c("estimate",
                                                       "std.error")])
  show("avg_comparisons x 1 to 2", avg_comparisons(
    f2, variables = list(x = c(1, 2))))
  show("emmeans at x = 1:2", emmeans(f2, ~ x, at = list(x = 1:2)))
  show("simo prior dirichlet(c(1, 1))", frm(
    y ~ mo(x) + z, data = d,
    prior = set_prior("dirichlet(c(1, 1))", class = "simo", coef = "mox1")))
  fp <- show("simo prior dirichlet(2)", frm(
    y ~ mo(x) + z, data = d,
    prior = set_prior("dirichlet(2)", class = "simo", coef = "mox1")))
  if (!inherits(fp, "error")) {
    cat("logLik with dirichlet(2):", format(as.numeric(logLik(fp)),
                                            digits = 12), "\n")
  }
  show("vcov names", rownames(vcov(f2)))
  show("confint", confint(f2))
  show("residuals head", head(residuals(f2)))
  show("simulate dim", dim(simulate(f2, nsim = 2, seed = 1)))
  show("update(f2) logLik", as.numeric(logLik(update(f2))))
}
fo <- show("fit y ~ mo(xo) + z (ordered factor, 2 levels)",
           frm(y ~ mo(xo) + z, data = d))
if (!inherits(fo, "error")) {
  cat("ordered logLik:", format(as.numeric(logLik(fo)), digits = 12), "\n")
  show("fixef ordered", fixef(fo))
}
fi <- show("fit y ~ mo(x) * z", frm(y ~ mo(x) * z, data = d))
if (!inherits(fi, "error")) {
  ml2 <- lm(y ~ x01 * z, data = d)
  cat("mo(x) * z logLik:", format(as.numeric(logLik(fi)), digits = 12),
      " lm:", format(as.numeric(logLik(ml2)), digits = 12), "\n")
  show("fixef mo(x) * z", fixef(fi))
}
fg <- show("fit y ~ mo(x) + z + (mo(x) | g)",
           frm(y ~ mo(x) + z + (mo(x) | g), data = d))
if (!inherits(fg, "error")) {
  cat("(mo(x) | g) logLik:", format(as.numeric(logLik(fg)), digits = 12),
      "\n")
  show("fixef (mo(x) | g)", fixef(fg))
  show("VarCorr (mo(x) | g)", VarCorr(fg))
}

## one value
f1 <- show("fit y ~ mo(x1) + z (one value)", frm(y ~ mo(x1) + z, data = d))
if (!inherits(f1, "error")) {
  ml1 <- lm(y ~ z, data = d)
  cat("one value logLik:", format(as.numeric(logLik(f1)), digits = 12),
      " lm(y ~ z):", format(as.numeric(logLik(ml1)), digits = 12), "\n")
  show("fixef one value", fixef(f1))
  show("summary one value", summary(f1))
  show("variables one value", variables(f1))
  show("default_prior one value", default_prior(y ~ mo(x1) + z, data = d))
  show("fitted newdata x1 = 4", fitted(f1, newdata = data.frame(x1 = 4,
                                                                z = 0)))
  show("fitted newdata x1 = 3", fitted(f1, newdata = data.frame(x1 = 3,
                                                                z = 0)))
  show("predictions one value", as.data.frame(predictions(
    f1, newdata = data.frame(x1 = 4, z = c(0, 1))))[, c("estimate",
                                                        "std.error")])
  show("conditional_effects z one value",
       conditional_effects(f1, effects = "z")$z[1:2, c("z", "x1",
                                                       "estimate__",
                                                       "se__")])
}
