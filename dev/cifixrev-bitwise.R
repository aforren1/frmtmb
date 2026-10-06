# Reviewer: outputs the cifix patch claims to leave unchanged to the
# bit, on fits that lost no standard error and on grids without a
# repeated gp() position. Reference BLAS.
# Usage: Rscript dev/cifixrev-bitwise.R <base|patch>   (writes
# dev/cifixrev-log/bitwise-<arm>.rds), then
# Rscript dev/cifixrev-bitwise.R compare
args <- commandArgs(TRUE)
arm <- args[1]
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-cifix/dev/cifixrev-log"
if (arm == "compare") {
  a <- readRDS(file.path(out_dir, "bitwise-base.rds"))
  b <- readRDS(file.path(out_dir, "bitwise-patch.rds"))
  n <- 0L
  same <- 0L
  for (m in names(a)) {
    for (o in names(a[[m]])) {
      n <- n + 1L
      x <- a[[m]][[o]]
      y <- b[[m]][[o]]
      ok <- identical(x, y)
      same <- same + ok
      if (!ok) {
        num <- function(z) {
          if (inherits(z, "Matrix")) z <- as.matrix(z)
          if (is.list(z)) z <- unlist(lapply(z, function(e) {
            if (is.numeric(e)) e else NULL
          }))
          as.numeric(z)
        }
        xx <- tryCatch(num(x), error = function(e) NULL)
        yy <- tryCatch(num(y), error = function(e) NULL)
        msg <- if (length(xx) && length(xx) == length(yy)) {
          k <- is.finite(xx) & is.finite(yy)
          sprintf("max |diff| / max |x| = %.3g, %d of %d entries differ",
                  max(abs(xx[k] - yy[k])) / max(abs(xx[k])),
                  sum(xx[k] != yy[k]), sum(k))
        } else "not numeric-comparable"
        cat("DIFFERS", m, o, ":", msg, "\n")
      }
    }
  }
  cat(sprintf("BITWISE %d of %d outputs identical\n", same, n))
  quit(save = "no")
}
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
libs <- c(base, user)
if (arm == "patch") libs <- c("C:/Users/adf44/source/r/cifixrev-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("arm", arm, find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
ns <- asNamespace("frmtmb")
q <- function(expr) suppressWarnings(expr)
out <- list()

## ordinal and mixture fits without random effects (dev/ordmix-bitwise.R)
set.seed(20261050)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
d$yn <- rnorm(n, ifelse(runif(n) < 0.4, 2, -1) + 0.5 * d$x)
om <- list(
  cum = function() frm(bf(y ~ x), family = cumulative(), data = d),
  cum_disc = function() frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(),
                            data = d),
  sratio_cs = function() frm(bf(y ~ cs(x)), family = sratio(), data = d),
  gauss_mix = function() frm(bf(yn ~ x), family = mixture(gaussian(),
                                                          gaussian()),
                             data = d))
for (m in names(om)) {
  fit <- q(om[[m]]())
  out[[m]] <- list(
    par = fit$opt$par, fitted = q(fitted(fit)),
    sim = simulate(fit, nsim = 2, seed = 1), fixef = q(fixef(fit)),
    vcov = q(vcov(fit)), jc = q(frm_joint_cov(fit)),
    linpred = q(frm_linpred(fit, newdata = d[1:5, ], se.fit = TRUE)),
    print = capture.output(print(fit)))
}

## random-effect fits
re_out <- function(fit, nd, ce = NULL) {
  r <- list(
    lost = ns$sdr_of(fit)$se_lost,
    par = fit$opt$par, vcov = q(vcov(fit)),
    jc = q(frm_joint_cov(fit)),
    lp_null = q(frm_linpred(fit, newdata = nd, se.fit = TRUE)),
    lp_na = q(frm_linpred(fit, newdata = nd, se.fit = TRUE,
                          re_formula = NA)),
    lb = q(frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)),
    ranef = q(ranef(fit, condVar = TRUE)), varcorr = q(VarCorr(fit)),
    fitted = q(fitted(fit, newdata = nd)),
    predict = {set.seed(3); q(predict(fit, newdata = nd))},
    summary = capture.output(print(q(summary(fit)))))
  if (!is.null(ce)) r$ce <- q(conditional_effects(fit, ce))
  r
}
ss <- lme4::sleepstudy
fs <- q(frm(Reaction ~ Days + (Days | Subject), data = ss))
out$sleep <- re_out(fs, ss[c(1, 15, 33), ], ce = "Days")
cb <- lme4::cbpp
fc <- q(frm(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial(), data = cb))
out$cbpp <- re_out(fc, cb[1:6, ], ce = "period")
out$cbpp$emm <- as.data.frame(q(emmeans::emmeans(fc, ~ period)))

set.seed(3)
d3 <- data.frame(x = stats::runif(200), z = stats::runif(200))
d3$y <- 1 + sin(4 * d3$x) + sin(2 * pi * d3$z) + stats::rnorm(200, 0, 0.3)
f3 <- q(frm(bf(y ~ s(x) + s(z)), data = d3))
out$smooth <- re_out(f3, data.frame(x = seq(0.05, 0.95, length.out = 9),
                                    z = 0.5), ce = "x")

## exact gp() at unseen positions, no position repeated
set.seed(17)
xg <- round(runif(120, 0, 10), 1)
dg <- data.frame(y = sin(xg) + rnorm(120, 0, 0.3), x = xg,
                 f = factor(sample(c("a", "b"), 120, TRUE)),
                 w = runif(120, 0.5, 2), x2 = runif(120, 0, 5))
dg$y2 <- dg$y + cos(dg$x2)
ndg <- data.frame(x = seq(0.05, 9.95, length.out = 41), f = "a",
                  w = 1.3, x2 = seq(0.1, 4.9, length.out = 41))
fg <- q(frm(bf(y ~ gp(x)), data = dg))
out$gp <- re_out(fg, ndg, ce = "x")
out$gp$curve <- q(frm_curve(fg, newdata = ndg[1:20, ], seed = 1))
fgb <- q(frm(bf(y ~ f + gp(x, by = f)), data = dg))
ndb <- rbind(ndg, transform(ndg, f = "b"))
out$gp_by <- re_out(fgb, ndb)
out$gp_by$diffcurve <- q(frm_curve(fgb, newdata = ndb[1:41, ],
                                   contrast = ndb[42:82, ],
                                   simultaneous = FALSE))
fg2 <- q(frm(bf(y2 ~ gp(x, x2)), data = dg))
out$gp2d <- re_out(fg2, ndg)
fgw <- q(frm(bf(y ~ gp(x, by = w)), data = dg))
out$gp_bynum <- re_out(fgw, transform(ndg, w = seq(0.6, 1.9,
                                                   length.out = 41)))
saveRDS(out, file.path(out_dir, paste0("bitwise-", arm, ".rds")))
cat("lost per fixture:\n")
for (m in names(out)) {
  cat(" ", m, ":", paste(names(out[[m]]$lost), collapse = ","), "\n")
}
