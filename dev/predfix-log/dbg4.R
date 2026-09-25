source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
set.seed(10)
db <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
db$y <- rbinom(60, 1, plogis(db$x + rnorm(10)[db$g]))
fq <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = db, quadrature = TRUE))
nd <- db[1:2, c("x", "g")]
r <- withCallingHandlers(frm_linpred(fq, newdata = nd, se.fit = TRUE, type = "response"),
  warning = function(w) { cat("W:", class(w)[1], conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(r$se.fit)
r <- withCallingHandlers(frm_linpred(fq, se.fit = TRUE)$se.fit[1:2],
  warning = function(w) { cat("W in-sample:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(r)
r <- withCallingHandlers(fitted(fq)[1:2,],
  warning = function(w) { cat("W fitted in-sample:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(r)
set.seed(1)
r <- withCallingHandlers(predict(fq, newdata = nd, ndraws=200),
  warning = function(w) { cat("W predict:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(r)
