# Review of lane loglik, attacks 2, 3 and 4: where logLik() refuses under
# REML, what AIC()/BIC() do for each kind of prior and each way of
# calling them, and what reaches AIC with no refusal.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
try_v <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste0("ERROR <", class(e)[1], "> ",
                       gsub("\\s+", " ", conditionMessage(e))),
                class = "err")
    }),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") },
    message = function(x) invokeRestart("muffleMessage"))
  if (inherits(v, "err")) return(unclass(v))
  paste0(paste(format(unclass(v), digits = 10), collapse = " "),
         if (length(w)) paste0("  [warnings: ", length(w), "]"))
}
say <- function(lbl, expr) cat(sprintf("%-46s %s\n", lbl, try_v(expr)))

set.seed(201)
d <- data.frame(x = rnorm(100), z = rnorm(100),
                g = factor(rep(1:20, each = 5)))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(100)

cat("\n## 2. REML refusal by class\n")
reml <- list(
  "sd only" = set_prior("normal(0, 1)", class = "sd"),
  "sigma only" = set_prior("exponential(1)", class = "sigma"),
  "b" = set_prior("normal(0, 1)", class = "b"),
  "Intercept, centered" = set_prior("normal(1, 1)", class = "Intercept"))
for (nm in names(reml)) {
  f <- suppressWarnings(frm(y ~ x + (1 | g), data = d, REML = TRUE,
                            prior = reml[[nm]]))
  say(paste("REML", nm, "logLik()"), as.numeric(logLik(f)))
  say(paste("REML", nm, "log_prior attr"),
      attr(logLik_penalized(f), "log_prior"))
}
f <- try(suppressWarnings(frm(bf(y ~ x + (1 | g), center = FALSE), data = d,
                              REML = TRUE,
                              prior = set_prior("normal(1, 1)", class = "b",
                                                coef = "Intercept"))))
if (!inherits(f, "try-error")) {
  say("REML Intercept center=FALSE (class b) logLik()", as.numeric(logLik(f)))
}
# distributional: sigma's coefficients under REML
f <- try(suppressWarnings(frm(bf(y ~ x + (1 | g), sigma ~ z), data = d,
                              REML = TRUE,
                              prior = set_prior("normal(0, 0.3)", class = "b",
                                                dpar = "sigma"))))
if (!inherits(f, "try-error")) {
  say("REML b dpar=sigma logLik()", as.numeric(logLik(f)))
  cat("   integrated comps:", frmtmb:::fit_integrated_comps(f), "\n")
}
# nonlinear under REML with a prior on an nlpar
set.seed(202)
dn <- data.frame(x = runif(80, 0, 4), g = factor(rep(1:10, each = 8)))
dn$y <- (3 + rnorm(10, 0, 0.3)[dn$g]) * exp(-0.7 * dn$x) + rnorm(80, 0, 0.15)
f <- try(suppressWarnings(frm(bf(y ~ a * exp(-k * x), a ~ 1 + (1 | g),
                                 k ~ 1, nl = TRUE), data = dn, REML = TRUE,
                              prior = set_prior("normal(1, 0.2)", nlpar = "k") +
                                set_prior("normal(3, 1)", nlpar = "a"))))
if (!inherits(f, "try-error")) {
  say("REML nlpar logLik()", as.numeric(logLik(f)))
} else cat("REML nlpar fit:", conditionMessage(attr(f, "condition")), "\n")
f <- suppressWarnings(frm(y ~ x + (1 | g), data = d, REML = TRUE,
                          prior = set_prior("normal(0, 1)", class = "b")))
say("REML b deviance()", deviance(f))
say("REML b summary()$loglik", summary(f)$loglik)
cat(grep("Method", capture.output(print(summary(f))), value = TRUE), "\n")

cat("\n## 3. AIC/BIC by kind of prior\n")
f0 <- frm(y ~ x + (1 | g), data = d)
kinds <- list(
  "none" = NULL,
  "explicit flat b ('')" = set_prior("", class = "b"),
  "bounds only b" = set_prior("", class = "b", coef = "x", lb = -5, ub = 5),
  "uniform(-5, 5) b" = set_prior("uniform(-5, 5)", class = "b", coef = "x"),
  "uniform(-5,5) b with lb/ub" = set_prior("uniform(-5, 5)", class = "b",
                                           coef = "x", lb = -5, ub = 5),
  "sd only gamma(2, 0.1)" = set_prior("gamma(2, 0.1)", class = "sd"),
  "normal(0, 1e6) b" = set_prior("normal(0, 1e6)", class = "b"))
for (nm in names(kinds)) {
  f <- try(suppressWarnings(frm(y ~ x + (1 | g), data = d,
                                prior = kinds[[nm]])), silent = TRUE)
  if (inherits(f, "try-error")) {
    cat(sprintf("%-46s FIT ERROR %s\n", nm,
                gsub("\\s+", " ", conditionMessage(attr(f, "condition")))))
    next
  }
  cat(sprintf("%-46s b_x %.8f  logLik %s  AIC %s\n", nm,
              f$estimates$beta[["x"]], try_v(as.numeric(logLik(f))),
              substr(try_v(AIC(f)), 1, 60)))
}
for (s in c("constant(0.5)", "flat", "horseshoe(1)", "R2D2()")) {
  say(paste0("prior string '", s, "'"),
      frm(y ~ x + (1 | g), data = d, prior = set_prior(s, class = "b")))
}

cat("\n## 3b. mixed calls\n")
fp <- frm(y ~ x + (1 | g), data = d,
          prior = set_prior("normal(0, 0.1)", class = "b"))
lmf <- lm(y ~ x, data = d)
say("AIC(fp)", AIC(fp))
say("AIC(f0, fp)", AIC(f0, fp)$AIC)
say("AIC(fp, k = log(100))", AIC(fp, k = log(100)))
say("AIC(lm, fp) (stats default first)", AIC(lmf, fp)$AIC)
say("BIC(lm, fp)", BIC(lmf, fp)$BIC)
say("AIC(logLik(fp))", AIC(logLik(fp)))
say("BIC(logLik(fp))", BIC(logLik(fp)))
say("-2 logLik(fp) + 2 df (what the two above use)",
    -2 * as.numeric(logLik(fp)) + 2 * attr(logLik(fp), "df"))
say("AIC(f0)", AIC(f0))
say("anova(f0, fp)", anova(f0, fp))
say("drop1(fp)", drop1(fp))
say("deviance(fp)", deviance(fp))
say("-2 * logLik(fp)", -2 * as.numeric(logLik(fp)))
fp2 <- frm(y ~ x + z + (1 | g), data = d,
           prior = set_prior("normal(0, 0.1)", class = "b"))
say("anova(fp, fp2)", anova(fp, fp2))
if (requireNamespace("performance", quietly = TRUE)) {
  say("performance::performance_aic(fp)", performance::performance_aic(fp))
}

cat("\n## 4. frm_multiple\n")
set.seed(203)
imps <- lapply(1:2, function(i) {
  dd <- data.frame(x = rnorm(50)); dd$y <- 1 + dd$x + rnorm(50); dd })
m1 <- frm_multiple(bf(y ~ x), data = imps,
                   prior = set_prior("normal(0, 1)", class = "b"))
say("logLik_penalized(frm_multiple)", logLik_penalized(m1))
say("logLik(frm_multiple)", logLik(m1))
say("logLik_penalized(lm)", logLik_penalized(lmf))
cat("\n## print of logLik_penalized\n")
print(logLik_penalized(fp))
print(logLik(fp))
cat("done\n")
