# A row of new data with NA in x or sdx, in every reader, beside a plain predictor
# with NA in z. Each answer, how many warnings the call raised, of which class.
#   Rscript dev/predread-me-na.R <library>
# Data seed 7421.
args <- commandArgs(TRUE)
.libPaths(c(args[1L], "C:/Users/adf44/source/r/rellib-r16", .libPaths()))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(7421)
n <- 200
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
d$y <- 1 + d$lat + 0.5 * d$z + rnorm(n, 0, 0.5)
d$cnt <- rpois(n, exp(0.3 + 0.6 * d$lat))
d$g <- factor(rep(1:20, each = 10))
d$xg <- as.numeric(tapply(d$lat, d$g, mean))[d$g]
d$sg <- 0.3
d$yg <- 1 + d$xg + rnorm(n, 0, 0.5)
fit <- suppressWarnings(frm(y ~ me(xo, sx) + z, data = d))
fp <- suppressWarnings(frm(cnt ~ me(xo, sx) + z, family = poisson(), data = d))
fg <- suppressWarnings(frm(yg ~ me(xg, sg, gr = g), data = d))
f0 <- suppressWarnings(frm(y ~ xo + z, data = d))
fz <- data.frame(xo = c(0, NA, 1, 1), sx = c(0.4, 0.4, NA, 0.4), z = c(0, 0, 0, NA))
show <- function(lab, e) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(e, warning = function(c) {
    w <<- c(w, class(c)[1L])
    invokeRestart("muffleWarning")
  }), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 70)))
  v <- if (is.character(r)) r else if (is.list(r)) {
    paste(signif(as.numeric(unlist(r[[1L]])), 3), collapse = " ")
  } else if (is.data.frame(r)) {
    paste(signif(as.numeric(r[[1L]]), 3), collapse = " ")
  } else paste(signif(as.numeric(r), 3), collapse = " ")
  cat(sprintf("%-28s %s\n%-28s warnings %d: %s\n", lab, substr(v, 1, 100), "",
              length(w), paste(unique(w), collapse = ",")))
}
run <- function(fit, nd, tag) {
  cat("---", tag, "\n")
  show("frm_linpred", frm_linpred(fit, newdata = nd))
  show("frm_linpred response", frm_linpred(fit, newdata = nd, type = "response"))
  show("frm_linpred se", frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit)
  show("frm_lp_basis extra", frm_lp_basis(fit, newdata = nd)$extra_var)
  show("fitted", fitted(fit, newdata = nd))
  set.seed(1); show("predict", predict(fit, newdata = nd, ndraws = 50))
  set.seed(1); show("simulate", simulate(fit, newdata = nd, nsim = 3))
  show("frm_curve", frm_curve(fit, newdata = nd, simultaneous = FALSE)$.estimate)
  show("emmeans", as.data.frame(emmeans::emmeans(fit, ~ z, at = list(z = 0))))
}
run(f0, transform(fz, xo = c(0, NA, 1, 1), sx = NULL), "plain, NA in xo and z")
run(fit, fz, "me gaussian")
run(fp, fz, "me poisson")
ndg <- data.frame(xg = c(0.1, 0.1, NA), sg = c(0.3, NA, 0.3), g = factor(c(1, 1, 2), 1:20))
run(fg, ndg, "me(gr = )")
cat("--- draws\n")
ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 200, seed = 3, refresh = 0))
show("posterior_epred", colMeans(posterior_epred(ds, newdata = fz)))
show("posterior_linpred", colMeans(posterior_linpred(ds, newdata = fz)))
show("posterior_predict", colMeans(posterior_predict(ds, newdata = fz)))
cat("--- the fit's own data, rows dropped for NA in x\n")
dn <- d
dn$xo[c(5, 9)] <- NA
fn <- suppressWarnings(frm(y ~ me(xo, sx) + z, data = dn))
show("fitted own data", nrow(fitted(fn, newdata = dn)))
show("NA rows", which(is.na(suppressWarnings(fitted(fn, newdata = dn))[, 1L])))
show("predict own data", nrow(predict(fn, newdata = dn, ndraws = 20)))
cat("--- a default draw stream without NA is unchanged\n")
set.seed(9)
print(signif(simulate(fit, newdata = fz[1, ], nsim = 2), 8))
