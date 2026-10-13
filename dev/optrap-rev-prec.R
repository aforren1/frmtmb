# Reviewer of lane optrap, item 8: how far the estimates of the rewritten
# precision tests are from their external references, in the fit's own
# standard errors, on one build. Run once per build and compare.
#   Rscript dev/optrap-rev-prec.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
for (e in parse("tests/testthat/test-autocor.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, env)
  }
}
out <- list()
put <- function(case, what, x, ref, se) {
  d <- abs(unname(x) - unname(ref))
  out[[length(out) + 1L]] <<- data.frame(
    build = lab, case = case, what = what, max_abs = max(d),
    max_se_units = max(d / se))
}
# IQ (test-bcm-gaussian.R)
x <- matrix(c(90, 95, 100, 105, 110, 115, 150, 155, 160), 3, byrow = TRUE)
d <- data.frame(y = as.vector(t(x)), person = factor(rep(1:3, each = 3)))
f <- frm(y ~ 0 + person, data = d)
co <- summary(f)$coefficients
put("IQ", "mu", f$estimates$beta, c(95, 110, 155), co$mu[, "Std. Error"])
sig <- sqrt(mean((d$y - rep(c(95, 110, 155), each = 3))^2))
put("IQ", "log sigma", f$estimates$betad[1], log(sig),
    co$sigma[, "Std. Error"])
# sleepstudy (test-glmm-gaussian.R)
f <- frm(Reaction ~ Days + (Days | Subject), data = sleepstudy)
l4 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
se <- summary(f)$coefficients$mu[, "Std. Error"]
put("sleepstudy", "fixef vs lmer", f$estimates$beta, fixef(l4), se)
put("sleepstudy", "logLik vs lmer", as.numeric(logLik(f)),
    as.numeric(logLik(l4)), 1)
gt <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy,
                       REML = FALSE)
put("sleepstudy", "fixef vs glmmTMB", f$estimates$beta,
    glmmTMB::fixef(gt)$cond, se)
put("sleepstudy", "SE / glmmTMB SE - 1", se /
      summary(gt)$coefficients$cond[, "Std. Error"], c(1, 1), 1)
# autocor (test-autocor.R)
cases <- list(
  ar1 = list(d = env$ac_sim(),
             fs = y ~ x + ar(week, subj, cov = TRUE),
             cor = function() nlme::corAR1(form = ~ week | subj)),
  ma1 = list(d = env$ac_sim(seed = 21, ar = numeric(0), ma = 0.7),
             fs = y ~ x + ma(week, subj, cov = TRUE),
             cor = function() nlme::corARMA(0.5, form = ~ week | subj,
                                            p = 0, q = 1)))
set.seed(5)
dc <- expand.grid(week = 1:5, subj = factor(1:40))
dc$x <- stats::rnorm(nrow(dc))
dc$y <- 1 + 0.5 * dc$x + stats::rnorm(40, 0, 1.2)[dc$subj] +
  stats::rnorm(nrow(dc))
cases$cosy <- list(d = dc, fs = y ~ x + cosy(week, subj),
                   cor = function() nlme::corCompSymm(form = ~ week | subj))
for (nm in names(cases)) for (m in c("ML", "REML")) {
  cs <- cases[[nm]]
  f <- frm(bf(cs$fs) + gaussian(), data = cs$d, REML = (m == "REML"))
  g <- nlme::gls(y ~ x, data = cs$d, correlation = cs$cor(), method = m,
                 control = env$gls_ctl())
  co <- summary(f)$coefficients
  k <- paste(nm, m)
  put(k, "logLik vs gls", as.numeric(logLik(f)), as.numeric(logLik(g)), 1)
  put(k, "beta", f$estimates$beta, coef(g), co$mu[, "Std. Error"])
  put(k, "log sigma", log(sigma(f)), log(g$sigma),
      co$sigma[, "Std. Error"])
  put(k, "cor", env$ac_pars(f),
      coef(g$modelStruct$corStruct, unconstrained = FALSE), env$ac_se(f))
}
X <- do.call(rbind, out)
write.table(X, sprintf("dev/optrap-rev-log/prec-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(X)
