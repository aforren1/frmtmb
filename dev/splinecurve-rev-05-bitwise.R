# Reviewer, lane splinecurve, claim 4: with allow_new_levels = FALSE every
# existing call is bitwise unchanged.
#   Rscript splinecurve-rev-05-bitwise.R <lib> <out.rds>
# Runs a set of curve calls modelled on the spline suite's fixtures and
# saves each result with its "fit" attribute removed (a fit holds a TMB
# pointer). The spec's new `allow_new_levels` field is kept here and
# removed only in the comparison, which says so.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
OUT <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline from", dirname(find.package("frmtmb.spline")), "\n")
out <- list()
keep <- function(name, expr) {
  ws <- character(0)
  v <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (!inherits(v, "err")) attr(v, "fit") <- NULL
  out[[name]] <<- list(value = v, warnings = ws)
  cat(name, if (inherits(v, "err")) "ERROR" else "ok", length(ws),
      "warnings\n")
}

## test-curve.R shapes
set.seed(1)
d1 <- data.frame(x = sort(runif(200)))
d1$y <- 2 * sin(pi * d1$x) + rnorm(200, 0, 0.4)
f1 <- frm(bf(y ~ s(x, k = 8)), family = gaussian(), data = d1)
g1 <- data.frame(x = seq(0.05, 0.95, length.out = 19))
keep("s_sim", frm_curve(f1, newdata = g1, seed = 1, nsim = 2000))
keep("s_deriv1", frm_curve_deriv(f1, var = "x", newdata = g1, seed = 2,
                                 nsim = 2000))
keep("s_deriv2", frm_curve_deriv(f1, var = "x", order = 2, newdata = g1,
                                 simultaneous = FALSE))
cv <- frm_curve(f1, newdata = g1, simultaneous = FALSE)
keep("s_stored_deriv", frm_curve_deriv(cv, var = "x", simultaneous = FALSE))
keep("s_feature_max", frm_curve_feature(f1, var = "x", type = "maximum",
                                        newdata = g1))
keep("s_feature_cross", frm_curve_feature(cv, var = "x", type = "crossing",
                                          at = 1))

set.seed(2)
dp <- data.frame(x = runif(300))
dp$y <- rpois(300, exp(1 + sin(2 * pi * dp$x)))
fp <- frm(bf(y ~ s(x, k = 8)), family = poisson(), data = dp)
keep("pois_transform", frm_curve(fp, newdata = g1, transform = TRUE,
                                 seed = 3, nsim = 2000))

set.seed(3)
dg <- data.frame(x = runif(300), g = factor(rep(1:15, each = 20)))
dg$y <- sin(2 * pi * dg$x) + rnorm(15, 0, 0.5)[dg$g] + rnorm(300, 0, 0.3)
fg <- frm(bf(y ~ s(x, k = 9) + (1 | g)), family = gaussian(), data = dg)
gg <- data.frame(x = seq(0, 1, length.out = 15),
                 g = factor("2", levels = levels(dg$g)))
keep("re_NA", frm_curve(fg, newdata = gg, seed = 4, nsim = 2000))
keep("re_NULL_seen", frm_curve(fg, newdata = gg, re_formula = NULL, seed = 4,
                               nsim = 2000))
keep("re_NULL_seen_contrast",
     frm_curve(fg, newdata = gg, re_formula = NULL, simultaneous = FALSE,
               contrast = transform(gg, g = factor("3",
                                                   levels = levels(dg$g)))))
keep("re_unseen_default",
     frm_curve(fg, newdata = transform(gg, g = factor("99",
                                                     levels = c(levels(dg$g),
                                                                "99"))),
               re_formula = NULL, simultaneous = FALSE))

set.seed(5)
ds <- data.frame(x = runif(300))
ds$y <- rnorm(300, sin(2 * pi * ds$x), exp(-1 + ds$x))
fsig <- frm(bf(y ~ x, sigma ~ s(x, k = 6)), family = gaussian(), data = ds)
keep("sigma_dpar", frm_curve(fsig, newdata = g1, dpar = "sigma",
                             simultaneous = FALSE))

## test-difference.R shapes
set.seed(6)
dd <- data.frame(x = runif(240), fac = factor(rep(c("A", "B"), 120)))
dd$y <- sin(2 * pi * dd$x) + (dd$fac == "B") * dd$x + rnorm(240, 0, 0.3)
fd <- frm(bf(y ~ fac + s(x, by = fac, k = 8)), family = gaussian(), data = dd)
gA <- data.frame(x = seq(0.05, 0.95, length.out = 12),
                 fac = factor("A", levels = c("A", "B")))
gB <- transform(gA, fac = factor("B", levels = c("A", "B")))
keep("diff", frm_curve(fd, newdata = gB, contrast = gA, seed = 7,
                       nsim = 2000))
cvd <- frm_curve(fd, newdata = gB, contrast = gA, simultaneous = FALSE)
keep("diff_deriv", frm_curve_deriv(cvd, var = "x", simultaneous = FALSE))
keep("diff_feature", frm_curve_feature(cvd, var = "x", type = "crossing",
                                       at = 0.3))

set.seed(8)
dgp <- data.frame(x = sort(runif(60, 0, 5)),
                  fac = factor(rep(c("A", "B"), 30)))
dgp$y <- sin(dgp$x) + ifelse(dgp$fac == "B", 0.5, 0) + rnorm(60, 0, 0.3)
fgp <- frm(bf(y ~ fac + gp(x)), family = gaussian(), data = dgp)
gx <- dgp$x[-1] - diff(dgp$x) / 2
gpA <- data.frame(x = gx, fac = factor("A", levels = levels(dgp$fac)))
gpB <- transform(gpA, fac = factor("B", levels = levels(dgp$fac)))
keep("gp_curve", frm_curve(fgp, newdata = gpA, simultaneous = FALSE))
keep("gp_diff", frm_curve(fgp, newdata = gpB, contrast = gpA,
                          simultaneous = FALSE))

## test-span.R shape: nl body with ps(), past the knot span
set.seed(4242)
dn <- data.frame(t = sort(runif(220)))
dn$y <- 2 + sin(2 * pi * dn$t) + rnorm(220, 0, 0.25)
fn <- frm(bf(y ~ lev + ps(t, k = 10, pad = 0.3), lev ~ 1, nl = TRUE),
          family = gaussian(), data = dn)
span <- fn$frame$linpreds[["y.mu"]]$ps_terms[[1]]$knot_range
keep("nl_inside", frm_curve(fn, newdata = data.frame(
  t = seq(span[1], span[2], length.out = 20)), seed = 9, nsim = 2000))
keep("nl_past_span", frm_curve(fn, newdata = data.frame(
  t = seq(span[2] - 0.2, span[2] + 1, length.out = 20)),
  simultaneous = FALSE))
keep("nl_deriv", frm_curve_deriv(fn, var = "t", newdata = data.frame(
  t = seq(span[1] + 0.05, span[2] - 0.05, length.out = 20)),
  simultaneous = FALSE))

## the fs model at a seen subject (test-curve.R, test-new-levels.R)
set.seed(4)
n_sub <- 12; n_rep <- 4; n_t <- 25
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
dfs <- data.frame(subject = factor(sub),
                  t = rep(seq(0, 1, length.out = n_t),
                          times = n_sub * n_rep))
h <- rnorm(n_sub, 1, 0.12)
sv <- rnorm(n_sub, 0.5, 0.04)
dfs$v <- peak(dfs$t, h[sub], sv[sub]) + rnorm(nrow(dfs), 0, 0.06)
ffs <- suppressWarnings(frm(bf(v ~ s(t, k = 10) +
                                 s(t, subject, bs = "fs", k = 5)),
                            family = gaussian(), data = dfs))
gfs <- data.frame(t = seq(0.05, 0.95, length.out = 30),
                  subject = factor("1", levels = levels(dfs$subject)))
keep("fs_seen", frm_curve(ffs, newdata = gfs, seed = 10, nsim = 2000))
keep("fs_seen_feature", frm_curve_feature(ffs, var = "t", type = "maximum",
                                          newdata = gfs))
saveRDS(out, OUT)
cat("saved", length(out), "results to", OUT, "\n")
