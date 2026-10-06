# Every post-fit method on ordinal mixtures, lane build. Seed 20261005.
# Usage: Rscript dev/ordmix-smoke.R [lane|base]
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb:", find.package("frmtmb"), "\n")
set.seed(20261005)
n <- 600
x <- rnorm(n)
z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.5 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
yh <- ifelse(runif(n) < 0.2, 0L, y)
d <- data.frame(y, yh, x, z, g)
nd <- d[1:5, ]
tryf <- function(label, expr, show = TRUE) {
  cat("\n=====", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
    if (show) print(r)
  invisible(r)
}
post <- function(f, d) {
  if (inherits(f, "error")) return(invisible())
  tryf("logLik", logLik(f))
  tryf("summary", summary(f))
  tryf("fixef", fixef(f))
  tryf("vcov dim", dim(vcov(f)))
  tryf("confint", confint(f))
  tryf("variables", variables(f))
  tryf("fitted head", head(fitted(f)[, 1, ]))
  tryf("fitted newdata", fitted(f, newdata = nd)[, 1, ])
  tryf("predict head", head(predict(f)))
  tryf("predict newdata", predict(f, newdata = nd))
  tryf("predict linear mu1", head(predict(f, dpar = "mu1", type = "link")))
  tryf("simulate", table(simulate(f, nsim = 1, seed = 1)[[1]]))
  tryf("simulate newdata", simulate(f, nsim = 2, seed = 1, newdata = nd))
  tryf("residuals", head(residuals(f)))
  tryf("residuals pearson", head(residuals(f, type = "pearson")))
  tryf("mixture_probs", head(mixture_probs(f)))
  tryf("log_lik", dim(log_lik(f)))
  tryf("log_lik sum vs logLik",
       c(sum(log_lik(f)), as.numeric(logLik(f))))
  tryf("conditional_effects x", {
    ce <- conditional_effects(f, effects = "x")
    head(ce[[1]])
  })
  tryf("conditional_effects categorical = FALSE", {
    ce <- conditional_effects(f, effects = "x", categorical = FALSE)
    head(ce[[1]])
  })
  tryf("hypothesis", hypothesis(f, "mu1_x > 0"))
  tryf("default_prior", default_prior(f))
  tryf("prior_summary", prior_summary(f))
  tryf("nobs", nobs(f))
  tryf("update", logLik(update(f, data = d)))
}
cat("\n############ order none, cumulative x2\n")
f1 <- tryf("fit", frm(bf(y ~ x), family = mixture(cumulative(),
                                                   cumulative()),
                      data = d), show = FALSE)
post(f1, d)
cat("\n############ order mu, cumulative x2\n")
f2 <- tryf("fit", frm(bf(y ~ x), family = mixture(cumulative(),
                                                   cumulative(),
                                                   order = "mu"),
                      data = d), show = FALSE)
post(f2, d)
cat("\n############ cumulative probit + sratio, disc1 ~ 0 + z, theta1 ~ z\n")
f3 <- tryf("fit", frm(bf(y ~ x, disc1 ~ 0 + z, theta1 ~ z),
                      family = mixture(cumulative("probit"), sratio()),
                      data = d), show = FALSE)
post(f3, d)
cat("\n############ thres(gr = g), sratio + acat\n")
f4 <- tryf("fit", frm(bf(y | thres(gr = g) ~ x),
                      family = mixture(sratio(), acat()),
                      data = d), show = FALSE)
post(f4, d)
cat("\n############ cs(x), sratio x2\n")
f5 <- tryf("fit", frm(bf(y ~ cs(x)), family = mixture(sratio(), sratio()),
                      data = d), show = FALSE)
post(f5, d)
cat("\n############ equidistant + sum_to_zero\n")
f6 <- tryf("fit", frm(bf(y ~ x),
                      family = mixture(cumulative(threshold = "equidistant"),
                                       cratio(threshold = "sum_to_zero")),
                      data = d), show = FALSE)
post(f6, d)
cat("\n############ hurdle x2\n")
f7 <- tryf("fit", frm(bf(yh ~ x), family = mixture(hurdle_cumulative(),
                                                    hurdle_cumulative()),
                      data = d), show = FALSE)
post(f7, d)
cat("\n############ (1 | g)\n")
f8 <- tryf("fit", frm(bf(y ~ x + (1 | g)),
                      family = mixture(cumulative(), cumulative()),
                      data = d), show = FALSE)
post(f8, d)
cat("\n############ order mu with priors, prior placement\n")
f9 <- tryf("fit", frm(bf(y ~ x), family = mixture(cumulative(), sratio(),
                                                   order = "mu"),
                      data = d,
                      prior = set_prior("normal(0, 1)", class = "Intercept")),
           show = FALSE)
post(f9, d)
cat("\n############ order none with prior on mu2 thresholds\n")
f10 <- tryf("fit", frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative()),
                       data = d,
                       prior = set_prior("normal(0, 1)", class = "Intercept",
                                         dpar = "mu2")), show = FALSE)
post(f10, d)
cat("\n############ refusals\n")
tryf("order mu + equidistant", frm(bf(y ~ x),
     family = mixture(cumulative(threshold = "equidistant"), cumulative(),
                      order = "mu"), data = d))
tryf("hurdle + cumulative", frm(bf(yh ~ x),
     family = mixture(hurdle_cumulative(), cumulative()), data = d))
tryf("groups =", frm(bf(y ~ x), family = mixture(cumulative(),
                                                  cumulative(), groups = ~g),
                     data = d))
tryf("ordinal + gaussian", mixture(cumulative(), gaussian()))
tryf("cs on cumulative component", frm(bf(y ~ cs(x)),
     family = mixture(cumulative(), sratio()), data = d))
tryf("cs on disc", frm(bf(y ~ x, disc ~ cs(z)), family = sratio(), data = d))
tryf("thres(gr) + cs", frm(bf(y | thres(gr = g) ~ cs(x)),
     family = mixture(sratio(), sratio()), data = d))
tryf("prior Intercept no dpar on order none",
     frm(bf(y ~ x), family = mixture(cumulative(), cumulative()), data = d,
         prior = set_prior("normal(0, 1)", class = "Intercept")))
tryf("prior Intercept dpar mu1 on order mu",
     frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                     order = "mu"), data = d,
         prior = set_prior("normal(0, 1)", class = "Intercept",
                           dpar = "mu1")))
cat("\n############ multivariate\n")
f11 <- tryf("fit", frm(bf(y ~ x) + bf(yh ~ x),
                       family = list(mixture(cumulative(), cumulative()),
                                     hurdle_cumulative()),
                       data = d), show = FALSE)
if (!inherits(f11, "error")) {
  tryf("summary", summary(f11))
  tryf("variables", variables(f11))
  tryf("fitted y", dim(fitted(f11, resp = "y")))
}
