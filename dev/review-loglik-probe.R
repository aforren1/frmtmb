# Review of lane loglik: probe the accessors the independent checks use.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
set.seed(11)
d <- data.frame(x = rnorm(80), g = factor(rep(1:10, each = 8)))
d$y <- 1 + 0.5 * d$x + rnorm(10, 0, 0.8)[d$g] + rnorm(10, 0, 0.3)[d$g] * d$x +
  rnorm(80, 0, 0.6)
fp <- frm(y ~ x + (1 + x | g), data = d,
          prior = set_prior("lkj(2)", class = "cor") +
            set_prior("normal(0, 1)", class = "sd"))
str(fp$estimates)
print(fixef(fp))
print(VarCorr(fp))
str(unclass(VarCorr(fp)))
print(sigma(fp))
print(head(fitted(fp)))
str(fp$frame$par_template)
print(logLik(fp)); print(logLik_penalized(fp))
ll <- tryCatch(log_lik(fp), error = function(e) conditionMessage(e))
str(ll)
