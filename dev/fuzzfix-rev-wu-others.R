# Reviewer: glmmTMB and lme4 on the models check_weights_unbounded()
# refuses (data of dev/fuzzfix-rev-wu.R, seed 20261021).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages({library(glmmTMB); library(lme4)})
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
cat(sprintf("sum(w) - n = %.4f\n", sum(d$w) - n))
show <- function(lbl, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-34s ERROR %s\n", lbl, conditionMessage(r)))
  } else {
    disp <- tryCatch(sigma(r), error = function(e) NA)
    cat(sprintf("%-34s logLik %.6g sigma() %.4g warnings %d %s\n", lbl,
                as.numeric(logLik(r)), disp, length(ws),
                if (length(ws)) substr(ws[1], 1, 90) else ""))
  }
}
show("glmmTMB gaussian (1|obs) w",
     glmmTMB(y ~ x + (1 | obs), d, weights = w))
show("glmmTMB gaussian (1|obs) w, unweighted",
     glmmTMB(y ~ x + (1 | obs), d))
show("glmmTMB Gamma (1|obs) w",
     glmmTMB(yg ~ x + (1 | obs), d, family = Gamma("log"), weights = w))
show("glmmTMB Gamma (1|obs) w*0.7",
     glmmTMB(yg ~ x + (1 | obs), d, family = Gamma("log"),
             weights = w * 0.7))
show("lmer (1|obs) w (precision weights)",
     lmer(y ~ x + (1 | obs), d, weights = w,
          control = lmerControl(check.nobs.vs.nlev = "ignore",
                                check.nobs.vs.nRE = "ignore")))
