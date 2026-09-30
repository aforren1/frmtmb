# Reviewer: the new refusals, each with its guard-absent case. Seed 91.
# REVLIB="" for the base arm.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  cat(sprintf("[%s] %s\n", label, if (inherits(r, "err")) paste("ERROR:", r)
              else "ok"))
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, unique(w)))
  invisible(if (inherits(r, "err")) NULL else r)
}
bf <- frmtmb::bf
set.seed(91)
n <- 50
d <- data.frame(x = rnorm(n), nt = rpois(n, 5) + 3L, w = runif(n, 1, 2))
d$yb <- rbinom(n, d$nt, 0.4)
d$y <- rnorm(n)
d$t <- rep(1:5, 10); d$g <- factor(rep(1:10, each = 5))
k <- 12L
d$yk <- rbinom(n, k, 0.4)
cat("\n-- a constant of the environment, bare and as an expression --\n")
tr("trials(k), k = 12 in env", frm(bf(yk | trials(k) ~ x), data = d,
                                   family = binomial()))
tr("trials(k + 0)", frm(bf(yk | trials(k + 0) ~ x), data = d,
                        family = binomial()))
tr("trials(12)", frm(bf(yk | trials(12) ~ x), data = d, family = binomial()))
tr("brms trials(k)", brms::standata(yk | trials(k) ~ x, data = d,
                                    family = stats::binomial()))
kw <- 2
tr("weights(kw) bare scalar", frm(bf(y | weights(kw) ~ x), data = d))
tr("brms weights(kw)", brms::standata(y | weights(kw) ~ x, data = d))
cat("\n-- wrong-length refusal and its absent case --\n")
tr("weights(w[1:3])", frm(bf(y | weights(w[1:3]) ~ x), data = d))
tr("weights(w[seq_along(w)])", frm(bf(y | weights(w[seq_along(w)]) ~ x),
                                    data = d))
tr("weights(rev(w))", frm(bf(y | weights(rev(w)) ~ x), data = d))
cat("\n-- 0 + intercept with an all-ones column (guard absent) --\n")
d1 <- d; d1$intercept <- 1
tr("intercept column of ones", frm(bf(y ~ 0 + intercept + x), data = d1))
d1$intercept <- TRUE
tr("intercept column TRUE", frm(bf(y ~ 0 + intercept + x), data = d1))
tr("brms intercept TRUE", brms::standata(y ~ 0 + intercept + x, data = d1))
cat("\n-- acformula resp matching (guard absent) --\n")
tr("bf(y ~ x) + acformula(resp = 'y')", frm(bf(y ~ x) +
                                              acformula(~ ar(t, g), resp = "y"),
                                            data = d))
tr("bf(y ~ x) + acformula(resp = 'z')", bf(y ~ x) +
     acformula(~ ar(t, g), resp = "z"))
tr("mv acformula resp = 'q'", bf(y ~ x) + bf(yb ~ x) +
     acformula(~ ar(t, g), resp = "q"))
cat("\n-- ARMA newdata response refusal and its absent cases --\n")
d$ya <- rnorm(n)
fa <- frm(bf(ya ~ x + ar(t, g)), data = d)
nd <- d[1:10, ]
nd$ya <- as.character(nd$ya)
tr("character response", predict(fa, newdata = nd))
nd$ya <- NA
tr("all-NA logical response", predict(fa, newdata = nd, ndraws = 5))
nd$ya <- NA_integer_
tr("all-NA integer response", predict(fa, newdata = nd, ndraws = 5))
nd$ya <- d$ya[1:10]; nd$ya[3] <- NaN
tr("NaN response", predict(fa, newdata = nd, ndraws = 5))
nd$ya <- d$ya[1:10]; nd$ya[3] <- Inf
tr("Inf response", print(fitted(fa, newdata = nd)[3:5, 1]))
cat("\n-- bernoulli refusals and absent cases --\n")
d$yb2 <- rbinom(n, 1, 0.5)
tr("bernoulli 0/1 integer", frm(bf(yb2 ~ x), data = d, family = bernoulli()))
d$yb3 <- d$yb2 + 0.0
tr("bernoulli 0/1 double", frm(bf(yb3 ~ x), data = d, family = bernoulli()))
d$yb4 <- d$yb2 * 1e10
tr("bernoulli 0/1e10", frm(bf(yb4 ~ x), data = d, family = bernoulli()))
d$yb5 <- d$yb2; d$yb5[2] <- 2
tr("bernoulli three values", frm(bf(yb5 ~ x), data = d, family = bernoulli()))
d$yb6 <- d$yb2; d$yb6[4] <- NA
tr("bernoulli with NA", frm(bf(yb6 ~ x), data = d, family = bernoulli()))
tr("bernoulli mi()", frm(bf(yb6 | mi() ~ x), data = d, family = bernoulli()))
