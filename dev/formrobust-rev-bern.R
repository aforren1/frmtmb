# Reviewer: bernoulli() coding of two values (claim 5). Data seed 51.
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
  if (inherits(r, "err")) cat(sprintf("[%s] ERROR: %s\n", label, r))
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, unique(w)))
  invisible(if (inherits(r, "err")) NULL else r)
}
same <- function(label, a, b) cat(sprintf("[%s] identical: %s\n", label,
                                          identical(a, b)))
sdY <- function(y, d) {
  d$yy <- y
  as.numeric(suppressMessages(suppressWarnings(
    brms::standata(yy ~ x, data = d, family = brms::bernoulli())))$Y)
}
set.seed(51)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, each = 10)))
d$y01 <- rbinom(n, 1, plogis(0.4 * d$x))
f01 <- frm(bf(y01 ~ x + (1 | g)), data = d, family = bernoulli())

cat("\n== 1. effect coding -0.5 / 0.5 ==\n")
d$yec <- d$y01 - 0.5
fe <- tr("-0.5/0.5", frm(bf(yec ~ x + (1 | g)), data = d,
                         family = bernoulli()))
if (!is.null(fe)) same("-0.5/0.5 logLik vs 0/1", logLik(fe), logLik(f01))
cat("brms Y for -0.5/0.5 equals 0/1:", identical(sdY(d$yec, d),
                                                 as.numeric(d$y01)), "\n")
d$y3 <- ifelse(d$y01 == 1, 3, 1)
f3 <- tr("1/3", frm(bf(y3 ~ x + (1 | g)), data = d, family = bernoulli()))
if (!is.null(f3)) same("1/3 logLik", logLik(f3), logLik(f01))
d$ypr <- runif(n)
tr("many fractions", frm(bf(ypr ~ x), data = d, family = bernoulli()))
tr("brms many fractions", sdY(d$ypr, d))

cat("\n== 2. -1 / -2 and the readers ==\n")
d$yn <- ifelse(d$y01 == 1, -1, -2)
fn <- frm(bf(yn ~ x + (1 | g)), data = d, family = bernoulli())
same("logLik", logLik(fn), logLik(f01))
same("simulate(seed)", simulate(fn, nsim = 2, seed = 3),
     simulate(f01, nsim = 2, seed = 3))
tr("refit with 0/1 codes", same("refit 0/1",
     logLik(refit(fn, simulate(fn, seed = 4)[[1]])),
     logLik(refit(f01, simulate(f01, seed = 4)[[1]]))))
tr("refit with original values", {
  r <- refit(fn, d$yn)
  cat("  refit(original -1/-2) logLik", format(logLik(r)), " fit",
      format(logLik(fn)), "\n") })
tr("bootstrap", same("frm_bootstrap", frm_bootstrap(fn, nsim = 3, seed = 5),
                     frm_bootstrap(f01, nsim = 3, seed = 5)))
tr("influence", same("influence fixed", influence(fn, groups = "g")$fixed,
                     influence(f01, groups = "g")$fixed))
nd <- d[d$yn == -2, ][1:6, ]
nd01 <- d[d$yn == -2, ][1:6, ]
tr("residuals newdata only -2", same("residuals(newdata only -2)",
     residuals(fn, newdata = nd), residuals(f01, newdata = nd01)))
tr("simulate newdata", same("simulate(newdata)",
     simulate(fn, nsim = 2, seed = 6, newdata = nd),
     simulate(f01, nsim = 2, seed = 6, newdata = nd01)))
tr("predict newdata", {
  set.seed(7); a <- predict(fn, newdata = nd)
  set.seed(7); b <- predict(f01, newdata = nd01)
  same("predict(newdata)", a, b) })
tr("update newdata only -2", {
  u <- update(fn, newdata = d[d$yn == -2 | seq_len(n) <= 3, ])
  cat("  update(newdata) ok, nobs", nobs(u), "\n") })
tr("pp_check newdata", { pp_check(fn, newdata = nd, ndraws = 3)
  cat("  pp_check ok\n") })
tr("residuals newdata with a third value", {
  nd3 <- nd; nd3$yn[1] <- 5; residuals(fn, newdata = nd3) })
tr("conditional_effects", same("ce estimates",
     conditional_effects(fn, effects = "x")[[1]]$estimate__,
     conditional_effects(f01, effects = "x")[[1]]$estimate__))
tr("bayes_R2", same("bayes_R2", bayes_R2(fn), bayes_R2(f01)))
tr("anova", invisible(anova(fn, update(fn, . ~ . - x))))

cat("\n== 3. logical and factor responses ==\n")
d$yl <- d$y01 == 1
fl <- frm(bf(yl ~ x + (1 | g)), data = d, family = bernoulli())
same("logical logLik", logLik(fl), logLik(f01))
d$ytrue <- TRUE
ft <- tr("all TRUE", frm(bf(ytrue ~ x), data = d, family = bernoulli()))
if (!is.null(ft)) cat("  all-TRUE intercept", fixef(ft)[1, 1], "\n")
cat("brms Y for all TRUE:", unique(sdY(d$ytrue, d)), "\n")
d$yone <- 5
fo <- tr("all 5", frm(bf(yone ~ x), data = d, family = bernoulli()))
if (!is.null(fo)) cat("  all-5 intercept", fixef(fo)[1, 1], "\n")
cat("brms Y for all 5:", unique(sdY(d$yone, d)), "\n")
d$yf1 <- factor("no", levels = c("no", "yes"))
ff1 <- tr("factor, one level used", frm(bf(yf1 ~ x), data = d,
                                        family = bernoulli()))
if (!is.null(ff1)) cat("  factor one-level intercept", fixef(ff1)[1, 1],
                       "\n")
cat("brms Y for factor one used level:", unique(sdY(d$yf1, d)), "\n")
d$ych <- ifelse(d$y01 == 1, "yes", "no")
fc <- frm(bf(ych ~ x + (1 | g)), data = d, family = bernoulli())
same("character logLik", logLik(fc), logLik(f01))

cat("\n== 4. multivariate, two bernoulli responses ==\n")
d$z01 <- rbinom(n, 1, 0.4)
d$zn <- ifelse(d$z01 == 1, "b", "a")
fm1 <- tr("mv coded", frm(bf(yn ~ x, family = bernoulli()) +
                            bf(zn ~ x, family = bernoulli()), data = d))
fm0 <- tr("mv 0/1", frm(bf(y01 ~ x, family = bernoulli()) +
                          bf(z01 ~ x, family = bernoulli()), data = d))
if (!is.null(fm1) && !is.null(fm0)) {
  same("mv logLik", logLik(fm1), logLik(fm0))
  same("mv simulate", unname(as.matrix(simulate(fm1, nsim = 1, seed = 8))),
       unname(as.matrix(simulate(fm0, nsim = 1, seed = 8))))
}
cat("\n== 5. bernoulli in a mixture or with mi / weights ==\n")
d$w <- runif(n, 0.5, 2)
fw <- tr("weights", frm(bf(yn | weights(w) ~ x), data = d,
                        family = bernoulli()))
fw0 <- frm(bf(y01 | weights(w) ~ x), data = d, family = bernoulli())
if (!is.null(fw)) same("weights logLik", logLik(fw), logLik(fw0))
fx <- tr("mixture(bernoulli, bernoulli)",
         frm(bf(yn ~ x), data = d, family = mixture(bernoulli(), bernoulli())))
fx0 <- tr("mixture 0/1", frm(bf(y01 ~ x), data = d,
                             family = mixture(bernoulli(), bernoulli())))
if (!is.null(fx) && !is.null(fx0)) same("mixture logLik", logLik(fx),
                                       logLik(fx0))
