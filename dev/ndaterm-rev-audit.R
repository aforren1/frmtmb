# Reviewer of lane ndaterm: every core newdata path on a binomial
# trials() fit, plus the legitimately-absent cases.
# Usage: Rscript dev/ndaterm-rev-audit.R <arm>   (mine | r9)
# Log: dev/ndaterm-rev-out/audit-<arm>.txt (redirected by the caller)
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste("ERROR:", conditionMessage(e)), class = "err")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [", label, "] WARN:", unique(w), "\n")
  if (inherits(r, "err")) cat("  [", label, "]", substr(r, 1, 300), "\n")
  invisible(r)
}
ok <- function(r) !inherits(r, "err") && !is.null(r)
show <- function(label, x) {
  cat(sprintf("%-50s %s\n", label,
              paste(format(x, digits = 5), collapse = "  ")))
}
est <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    if ("Estimate" %in% colnames(x)) return(x[, "Estimate"])
    if ("estimate" %in% colnames(x)) return(x[, "estimate"])
    if ("Predicted" %in% colnames(x)) return(x[, "Predicted"])
    return(colMeans(as.matrix(x)))
  }
  as.numeric(x)
}

set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE),
                w = runif(n, 0.5, 2), g = gl(10, 20))
d$y <- rbinom(n, d$N, plogis(0.3 + 0.5 * d$x))
d$y10 <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x))
d$M <- d$N - 1
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
nd <- data.frame(x = c(0, 0), N = c(1000, 2), w = 1, g = "1", M = c(999, 1))

cat("\n== A. binomial trials(N): every core newdata path\n")
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = d)
f0 <- fitted(fb, newdata = nd)
show("fitted(newdata) [reference]", est(f0))
show("simulate(newdata) row means",
     rowMeans(as.matrix(simulate(fb, nsim = 400, newdata = nd, seed = 1))))
set.seed(1)
show("predict(newdata) Estimate", est(predict(fb, newdata = nd,
                                              ndraws = 400)))
r <- tr("frm_linpred response", frm_linpred(fb, newdata = nd,
                                             type = "response"))
if (ok(r)) show("frm_linpred(type = response)", est(r))
r <- tr("predict se.fit", predict(fb, newdata = nd, type = "response",
                                  se.fit = TRUE))
if (ok(r)) show("predict(se.fit = TRUE)", unlist(r)[1:2])
r <- tr("fitted re_formula NA", fitted(fb, newdata = nd, re_formula = NA))
if (ok(r)) show("fitted(re_formula = NA)", est(r))
ndy <- transform(nd, y = c(560, 1))
r <- tr("residuals(newdata)", residuals(fb, newdata = ndy))
if (ok(r)) show("residuals(newdata) (y - fitted)", est(r))
r <- tr("pp_check(newdata)", suppressMessages(
  pp_check(fb, newdata = ndy, ndraws = 20, type = "stat")))
if (ok(r)) cat("pp_check(newdata) ok\n")
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  r <- tr("marginaleffects::predictions",
          marginaleffects::predictions(fb, newdata = nd))
  if (ok(r)) show("marginaleffects predictions", r$estimate)
}
if (requireNamespace("insight", quietly = TRUE)) {
  r <- tr("insight::get_predicted", insight::get_predicted(fb, data = nd))
  if (ok(r)) show("insight get_predicted", as.numeric(r))
  r <- tr("insight::get_predicted prediction",
          insight::get_predicted(fb, data = nd, predict = "prediction"))
  if (ok(r)) show("insight get_predicted(prediction)", as.numeric(r))
}
if (requireNamespace("emmeans", quietly = TRUE)) {
  r <- tr("emmeans epred", summary(emmeans::emmeans(
    fb, ~ x, at = list(x = 0, N = 1000), epred = TRUE)))
  if (ok(r)) show("emmeans(epred, N = 1000)", r$emmean)
  r <- tr("emmeans epred default", summary(emmeans::emmeans(
    fb, ~ x, at = list(x = 0), epred = TRUE)))
  if (ok(r)) show("emmeans(epred, N default)", r$emmean)
}
set.seed(2)
r <- tr("ce epred N=1000", conditional_effects(
  fb, "x", conditions = data.frame(N = 1000))[[1]])
if (ok(r)) show("CE epred (N = 1000) at x~0",
                r$estimate__[which.min(abs(r$x))])
set.seed(2)
r <- tr("ce predict N=1000", conditional_effects(
  fb, "x", conditions = data.frame(N = 1000), method = "predict",
  ndraws = 201)[[1]])
if (ok(r)) show("CE predict (N = 1000) at x~0",
                r$estimate__[which.min(abs(r$x))])
set.seed(2)
r <- tr("ce epred default", conditional_effects(fb, "x")[[1]])
if (ok(r)) show("CE epred default trials at x~0 (brms: 1)",
                r$estimate__[which.min(abs(r$x))])
r <- tr("ce predict default", conditional_effects(fb, "x",
  method = "predict", ndraws = 201)[[1]])
if (ok(r)) show("CE predict default at x~0",
                r$estimate__[which.min(abs(r$x))])
r <- tr("dharma", dharma_residuals(fb, nsim = 50))
if (ok(r)) cat("dharma_residuals in-sample ok\n")

cat("\n== B. trials() constant and expression\n")
fc <- frm(bf(y10 | trials(10) ~ x), family = binomial(), data = d)
r <- tr("trials(10) predict nd without vars",
        predict(fc, newdata = data.frame(x = 0), ndraws = 200,
                summary = FALSE))
if (ok(r)) show("trials(10) predict col mean, max", c(mean(r), max(r)))
r <- tr("trials(10) fitted", fitted(fc, newdata = data.frame(x = 0)))
if (ok(r)) show("trials(10) fitted", est(r))
r <- tr("trials(10) simulate", simulate(fc, nsim = 200,
                                        newdata = data.frame(x = 0),
                                        seed = 3))
if (ok(r)) show("trials(10) simulate mean", mean(as.matrix(r)))
fe <- frm(bf(y | trials(M + 1) ~ x), family = binomial(), data = d)
set.seed(3)
r <- tr("trials(M+1) predict", predict(fe, newdata = nd, ndraws = 400,
                                       summary = FALSE))
if (ok(r)) show("trials(M + 1) predict col means", colMeans(r))
r <- tr("trials(M+1) fitted", fitted(fe, newdata = nd))
if (ok(r)) show("trials(M + 1) fitted", est(r))
r <- tr("trials(M+1) no M", predict(fe, newdata = data.frame(x = 0),
                                    ndraws = 5))
if (ok(r)) show("trials(M + 1) newdata without M", est(r))

cat("\n== C. weights(), cens(), mi(), subset()\n")
fw <- frm(bf(y | trials(N) + weights(w) ~ x), family = binomial(), data = d)
r <- tr("weights: nd without w", predict(fw, newdata = nd[, c("x", "N")],
                                         ndraws = 100))
if (ok(r)) show("weights fit, newdata without w", est(r))
r <- tr("weights: fitted nd without w",
        fitted(fw, newdata = nd[, c("x", "N")]))
if (ok(r)) show("weights fit fitted, no w", est(r))
r <- tr("weights: simulate nd without w",
        simulate(fw, nsim = 50, newdata = nd[, c("x", "N")], seed = 4))
if (ok(r)) show("weights fit simulate no w", rowMeans(as.matrix(r)))
dcn <- d
dcn$cc <- ifelse(dcn$cnt > 4, "right", "none")
dcn$cnt2 <- pmin(dcn$cnt, 4)
fcn <- frm(bf(cnt2 | cens(cc) ~ x), family = poisson(), data = dcn)
r <- tr("cens: predict nd without cc",
        predict(fcn, newdata = data.frame(x = c(0, 1)), ndraws = 100))
if (ok(r)) show("cens fit predict, no cc", est(r))
r <- tr("cens: predict nd with cc",
        predict(fcn, newdata = data.frame(x = c(0, 1), cc = "right"),
                ndraws = 100))
if (ok(r)) show("cens fit predict, with cc", est(r))
dmi <- d
dmi$gy <- 1 + dmi$x + rnorm(n)
dmi$gy[1:10] <- NA
fmi <- tr("mi fit", frm(bf(gy | mi() ~ x), family = gaussian(), data = dmi))
if (ok(fmi)) {
  r <- tr("mi predict nd", predict(fmi, newdata = data.frame(x = c(0, 1)),
                                   ndraws = 100))
  if (ok(r)) show("mi() fit predict(newdata)", est(r))
}
dsb <- d
dsb$sb <- rep(c(TRUE, FALSE), n / 2)
dsb$gy <- 1 + dsb$x + rnorm(n)
fsb <- tr("subset mv fit", frm(
  mvbf(bf(y | trials(N) + subset(sb) ~ x), bf(gy ~ x)),
  family = list(binomial(), gaussian()), data = dsb))
if (ok(fsb)) {
  r <- tr("subset mv predict nd with sb",
          predict(fsb, newdata = transform(nd, sb = TRUE), ndraws = 200))
  if (ok(r)) { cat("subset mv predict(newdata):\n"); print(r) }
  r <- tr("subset mv predict resp gy, nd without N or sb",
          predict(fsb, newdata = data.frame(x = 0), resp = "gy",
                  ndraws = 50))
  if (ok(r)) show("mv predict resp = gy, newdata no N", est(r))
}

cat("\n== D. legitimately absent\n")
r <- tr("in-sample predict", predict(fb, summary = FALSE, ndraws = 20))
if (ok(r)) cat("in-sample predict ok, max/N <= 1:",
               all(sweep(r, 2, d$N, "<=")), "\n")
r <- tr("in-sample predict re_formula NA",
        predict(fb, re_formula = NA, ndraws = 20))
if (ok(r)) cat("in-sample predict(re_formula = NA) ok\n")
fbg <- frm(bf(y | trials(N) ~ x + (1 | g)), family = binomial(), data = d)
r <- tr("re fit nd re_formula NA",
        predict(fbg, newdata = nd[, c("x", "N")], re_formula = NA,
                ndraws = 200))
if (ok(r)) show("RE fit predict(newdata, re_formula = NA)", est(r))
r <- tr("re fit nd", predict(fbg, newdata = nd, ndraws = 200))
if (ok(r)) show("RE fit predict(newdata)", est(r))
r <- tr("re fit in-sample re_formula NA",
        predict(fbg, re_formula = NA, ndraws = 20, summary = FALSE))
if (ok(r)) cat("RE fit in-sample re_formula = NA ok\n")
fp <- frm(bf(cnt ~ x), family = poisson(), data = d)
r <- tr("poisson nd with extra cols", predict(fp, newdata = nd,
                                              ndraws = 20))
if (ok(r)) cat("poisson predict(newdata) ok\n")
fpw <- frm(bf(cnt | weights(w) ~ x), family = poisson(), data = d)
r <- tr("poisson weights nd without w",
        predict(fpw, newdata = data.frame(x = 0), ndraws = 20))
if (ok(r)) cat("poisson weights predict(newdata without w) ok\n")
fz <- tr("zib fit", frm(bf(y | trials(N) ~ x), family =
                          zero_inflated_binomial(), data = d))
if (ok(fz)) {
  r <- tr("zib in-sample", predict(fz, ndraws = 20))
  if (ok(r)) cat("zero_inflated_binomial in-sample ok\n")
  r <- tr("zib nd", predict(fz, newdata = nd, ndraws = 400))
  if (ok(r)) show("zib predict(newdata) vs fitted",
                  c(est(r), est(fitted(fz, newdata = nd))))
}
fzbb <- tr("zibb fit", frm(bf(y | trials(N) ~ x), family =
                             zero_inflated_beta_binomial(), data = d))
if (ok(fzbb)) {
  r <- tr("zibb nd", predict(fzbb, newdata = nd, ndraws = 400))
  if (ok(r)) show("zibb predict(newdata) vs fitted",
                  c(est(r), est(fitted(fzbb, newdata = nd))))
}
fmx <- tr("mixture binomial fit", frm(
  bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
  family = mixture(binomial(), binomial()), data = d))
if (ok(fmx)) {
  r <- tr("mixture in-sample", predict(fmx, ndraws = 20))
  if (ok(r)) cat("mixture(binomial, binomial) in-sample ok\n")
  set.seed(5)
  r <- tr("mixture nd", predict(fmx, newdata = nd, ndraws = 400))
  if (ok(r)) show("mixture predict(newdata) vs fitted",
                  c(est(r), est(fitted(fmx, newdata = nd))))
  r <- tr("mixture nd without N", predict(fmx, newdata = data.frame(x = 0),
                                          ndraws = 5))
  if (ok(r)) show("mixture newdata without N", est(r))
}
fmxp <- tr("mixture binomial+poisson", frm(
  bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
  family = mixture(binomial(), poisson()), data = d))
if (ok(fmxp)) {
  r <- tr("mixture b+p nd", predict(fmxp, newdata = nd, ndraws = 400))
  if (ok(r)) show("mixture(binomial, poisson) predict vs fitted",
                  c(est(r), est(fitted(fmxp, newdata = nd))))
}
dmv <- d
dmv$gy <- 1 + dmv$x + rnorm(n)
fmv <- frm(mvbf(bf(y | trials(N) ~ x), bf(gy ~ x)),
           family = list(binomial(), gaussian()), data = dmv)
r <- tr("mv predict nd", predict(fmv, newdata = nd, ndraws = 400))
if (ok(r)) { cat("mvbf predict(newdata):\n"); print(r) }
r <- tr("mv predict resp gy nd without N",
        predict(fmv, newdata = data.frame(x = 0), resp = "gy", ndraws = 50))
if (ok(r)) show("mvbf predict(resp = gy, no N)", est(r))
r <- tr("mv fitted resp gy nd without N",
        fitted(fmv, newdata = data.frame(x = 0), resp = "gy"))
if (ok(r)) show("mvbf fitted(resp = gy, no N)", est(r))
r <- tr("mv simulate nd without N (all resp)",
        simulate(fmv, nsim = 5, newdata = data.frame(x = 0), seed = 1))
r <- tr("mv in-sample predict", predict(fmv, ndraws = 20))
if (ok(r)) cat("mvbf in-sample predict ok\n")

cat("\n== E. refusal texts\n")
r <- tr("binomial nd without N predict", predict(fb, newdata =
                                                   data.frame(x = 0)))
r <- tr("binomial nd without N fitted", fitted(fb, newdata =
                                                 data.frame(x = 0)))
r <- tr("binomial nd without N simulate", simulate(fb, newdata =
                                                     data.frame(x = 0)))
r <- tr("binomial nd without N frm_linpred response",
        frm_linpred(fb, newdata = data.frame(x = 0), type = "response"))
cat("done\n")
