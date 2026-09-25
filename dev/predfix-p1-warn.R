# Punch round 1, B1: does an engaged default fit warn more than
# autoscale = FALSE on the same data? The reviewer's constructions.
#   PREDFIX_ARM=lane Rscript dev/predfix-p1-warn.R > dev/predfix-log/p1-warn-<tag>.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
catch <- function(expr) {
  w <- character(0)
  f <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(fit = f, w = w)
}
row <- function(lbl, mk) {
  a <- catch(mk(NULL))
  b <- catch(mk(FALSE))
  pa <- a$fit$opt$par
  pb <- b$fit$opt$par
  cat(sprintf("%-44s default %d warn, FALSE %d warn; logLik %.9f / %.9f; engaged %s; max|dpar| %.3g\n",
              lbl, length(a$w), length(b$w), as.numeric(logLik(a$fit)),
              as.numeric(logLik(b$fit)), !is.null(a$fit$par_units),
              max(abs(pa - pb))))
  for (x in setdiff(a$w, b$w)) cat("    only under default:", substr(x, 1, 110), "\n")
}
ctl <- function(a) frmtmb_control(autoscale = a)

mk_skew <- function(seed, s, n = 250) {
  set.seed(seed)
  xs <- -abs(rnorm(n)) * 3
  yy <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = yy, xs = xs * s, xr = xs)
}
fo <- bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs)
for (sc in list(c(2, 5e-4), c(3, 5e-4), c(3, 9.9e-4), c(3, 1.01e-3),
                c(1, 1e-6))) {
  d <- mk_skew(sc[1], sc[2])
  row(sprintf("skew_normal alpha ~ 0 + xs, seed %d s %g", sc[1], sc[2]),
      function(a) frm(fo, family = skew_normal(), data = d, control = ctl(a)))
}

set.seed(1)
n <- 200
x <- rnorm(n) * 1e-4 + 5e-4
dg <- data.frame(x = x, y = 1 + 2000 * (x - 5e-4) + rnorm(n))
row("gaussian, Intercept ub = 0.5",
    function(a) frm(bf(y ~ x), family = gaussian(), data = dg,
                    prior = set_prior("", class = "Intercept", ub = 0.5),
                    control = ctl(a)))
row("gaussian, b(x) lb = 3e4",
    function(a) frm(bf(y ~ x), family = gaussian(), data = dg,
                    prior = set_prior("", class = "b", coef = "x", lb = 3e4),
                    control = ctl(a)))
