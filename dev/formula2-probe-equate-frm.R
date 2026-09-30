# The equated mixture through frmtmb's post-fit surface. Seed 11.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  cat("==", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1L && startsWith(r, "ERROR")) {
    cat(r, "\n")
  } else {
    print(r)
  }
  invisible(r)
}
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
fam <- mixture(gaussian(), gaussian())
fit <- frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)
tr("spec", fit$spec)
tr("logLik", logLik(fit))
tr("fixef", fixef(fit))
tr("variables", variables(fit))
tr("summary", summary(fit))
tr("vcov dim", dim(vcov(fit)))
tr("confint", confint(fit))
tr("get_prior", get_prior(bf(y ~ x, sigma1 = "sigma2"), family = fam,
                          data = d))
tr("prior on sigma1", frm(bf(y ~ x, sigma1 = "sigma2"), family = fam,
                          data = d, prior = set_prior("normal(1, 1)",
                                                      class = "sigma1")))
tr("prior on sigma2", logLik(frm(bf(y ~ x, sigma1 = "sigma2"),
                                 family = fam, data = d,
                                 prior = set_prior("normal(1, 1)",
                                                   class = "sigma2"))))
tr("fitted head", head(fitted(fit)))
tr("predict head", head(predict(fit)))
tr("predict newdata", predict(fit, newdata = d[1:3, ]))
tr("simulate", dim(as.matrix(simulate(fit, nsim = 2, seed = 1))))
tr("hypothesis", hypothesis(fit, "sigma1 = sigma2"))
tr("par_template", par_template(fit))
tr("frm_linpred sigma1 vs sigma2",
   all.equal(frm_linpred(fit, dpar = "sigma1"),
             frm_linpred(fit, dpar = "sigma2")))
tr("AIC df", attr(logLik(fit), "df"))
tr("update to lf spelling", logLik(update(fit, formula. = bf(y ~ x) +
                                            lf(sigma1 = "sigma2"))))
tr("print bf", bf(y ~ x, sigma1 = "sigma2"))
tr("equate to missing", frm(bf(y ~ x, sigma1 = "sigma3"), family = fam,
                            data = d))
tr("equate in gaussian", frm(bf(y ~ x, sigma = "sigma"), data = d))
tr("theta", bf(y ~ x, theta1 = "theta2"))
