# Reviewer of lane ndaterm: mixture row alignment of trials(), cens()
# on newdata, and dec() inputs on a wiener fit.
# Usage: Rscript dev/ndaterm-rev-probe2.R <arm>   (mine | r9)
# Log: dev/ndaterm-rev-out/probe2-<arm>.txt (redirected by the caller)
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "eam", find.package("frmtmb.eam"), "\n")
tr <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste("ERROR:", conditionMessage(e)), class = "err")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [", label, "] WARN:", substr(unique(w), 1, 150),
                     "\n")
  if (inherits(r, "err")) cat("  [", label, "]", substr(r, 1, 300), "\n")
  invisible(r)
}
ok <- function(r) !inherits(r, "err") && !is.null(r)
show <- function(label, x) {
  cat(sprintf("%-52s %s\n", label,
              paste(format(x, digits = 5), collapse = "  ")))
}

cat("\n== 1. mixture(binomial, binomial), per-row trials()\n")
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(c(2, 50), n, TRUE))
cls <- rbinom(n, 1, 0.5)
d$y <- rbinom(n, d$N, ifelse(cls == 1, plogis(-1 + 0.5 * d$x),
                             plogis(1 + 0.5 * d$x)))
fmx <- frm(bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
           family = mixture(binomial(), binomial()), data = d)
s <- tr("simulate in-sample", as.matrix(simulate(fmx, nsim = 200,
                                                 seed = 1)))
if (ok(s)) {
  show("simulate in-sample: share of draws > N", mean(s > d$N))
  show("simulate in-sample: mean draw at N = 2 rows",
       mean(s[d$N == 2, ]))
  show("fitted at N = 2 rows (mean)", mean(fitted(fmx)[d$N == 2,
                                                     "Estimate"]))
}
set.seed(2)
p <- tr("predict in-sample", predict(fmx, summary = FALSE, ndraws = 200))
if (ok(p)) show("predict in-sample: share of draws > N",
                mean(sweep(p, 2, d$N, ">")))
nd <- data.frame(x = c(0, 0), N = c(1000, 2))
set.seed(3)
p <- tr("predict newdata", predict(fmx, newdata = nd, summary = FALSE,
                                   ndraws = 400))
if (ok(p)) {
  show("predict(newdata) col means", colMeans(p))
  show("predict(newdata) col max (N = 1000, 2)", apply(p, 2, max))
}
s <- tr("simulate newdata", as.matrix(simulate(fmx, nsim = 400,
                                               newdata = nd, seed = 3)))
if (ok(s)) {
  show("simulate(newdata) row means", rowMeans(s))
  show("simulate(newdata) row max", apply(s, 1, max))
}
show("fitted(newdata)", fitted(fmx, newdata = nd)[, "Estimate"])

cat("\n== 2. poisson cens(cc): predict(newdata) with and without cc\n")
dc <- data.frame(x = rnorm(n))
dc$cnt <- rpois(n, exp(0.5 + 0.3 * dc$x))
dc$cc <- ifelse(dc$cnt > 4, "right", "none")
dc$cnt2 <- pmin(dc$cnt, 4)
fcn <- suppressMessages(frm(bf(cnt2 | cens(cc) ~ x), family = poisson(),
                            data = dc))
set.seed(4)
a <- predict(fcn, newdata = data.frame(x = c(0, 1)), ndraws = 4000)
set.seed(4)
b <- predict(fcn, newdata = data.frame(x = c(0, 1), cc = "right"),
             ndraws = 4000)
show("cens predict(newdata) without cc", a[, "Estimate"])
show("cens predict(newdata) with cc = right", b[, "Estimate"])
show("fitted(newdata)", fitted(fcn, newdata = data.frame(x = c(0, 1)))[,
                                                           "Estimate"])

cat("\n== 3. wiener dec() inputs\n")
set.seed(7)
dw <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.3)
dw$x <- rnorm(300)
dw$side <- factor(ifelse(dw$upper == 1, "up", "down"),
                  levels = c("up", "down"))
dw$lg <- dw$upper == 1
fw <- frm(bf(rt | dec(upper) ~ x, bias = 0.3), family = wiener(),
          data = dw)
ff <- frm(bf(rt | dec(side) ~ x, bias = 0.3), family = wiener(),
          data = dw)
fl <- tr("logical dec fit", frm(bf(rt | dec(lg) ~ x, bias = 0.3),
                                family = wiener(), data = dw))
ref <- fitted(fw, newdata = data.frame(x = 0, upper = c(0, 1)))[,
                                                         "Estimate"]
show("numeric fit fitted(upper = 0, 1) [reference]", ref)
r <- tr("CE default epred (numeric dec)", conditional_effects(fw, "x")[[1]])
if (ok(r)) show("CE default epred numeric dec at x~0",
                r$estimate__[which.min(abs(r$x))])
r <- tr("CE default (no effects arg)", conditional_effects(fw))
if (ok(r)) cat("CE default call ok\n")
r <- tr("CE epred upper = 1", conditional_effects(
  fw, "x", conditions = data.frame(upper = 1))[[1]])
if (ok(r)) show("CE epred conditions upper = 1",
                r$estimate__[which.min(abs(r$x))])
r <- tr("CE factor default", conditional_effects(ff, "x")[[1]])
if (ok(r)) show("CE epred factor default (first level up = 0)",
                r$estimate__[which.min(abs(r$x))])
r <- tr("CE over the dec factor itself", conditional_effects(ff, "side"))
if (ok(r)) show("CE over side: estimates", r[[1]]$estimate__)
r <- tr("CE over the numeric dec var itself", conditional_effects(fw,
                                                                  "upper"))
if (ok(r)) show("CE over upper: estimates", r[[1]]$estimate__)
fac <- function(v, lv = c("up", "down")) factor(v, levels = lv)
cases <- list(
  "factor both levels, fit order" = data.frame(x = 0, side = fac(c("down", "up"))),
  "factor one level present, both declared" = data.frame(x = 0, side = fac("down")),
  "factor dropped to one level" = data.frame(x = 0, side = factor("down")),
  "relevel(side, 'down')" = data.frame(x = 0, side = relevel(fac(c("down", "up")), "down")),
  "character c('down','up')" = data.frame(x = 0, side = c("down", "up")),
  "character 'down' alone" = data.frame(x = 0, side = "down"),
  "unseen levels left/right" = data.frame(x = 0, side = factor(c("left", "right"))),
  "factor with NA" = data.frame(x = 0, side = fac(c("down", NA))))
for (nm in names(cases)) {
  r <- tr(nm, fitted(ff, newdata = cases[[nm]]))
  if (ok(r)) show(paste("factor fit:", nm), r[, "Estimate"])
}
show("factor fit reference: down (= upper 1), up (= upper 0)",
     rev(ref))
if (ok(fl)) {
  for (v in list(c(FALSE, TRUE), TRUE, c(TRUE, NA))) {
    r <- tr(paste("logical", paste(v, collapse = ",")),
            fitted(fl, newdata = data.frame(x = 0, lg = v)))
    if (ok(r)) show(paste("logical fit:", paste(v, collapse = ",")),
                    r[, "Estimate"])
  }
}
for (v in list(c(0, 1), 2, -1, 0.5, c(1, NA), TRUE)) {
  r <- tr(paste("numeric fit upper =", paste(v, collapse = ",")),
          fitted(fw, newdata = data.frame(x = 0, upper = v)))
  if (ok(r)) show(paste("numeric fit upper =", paste(v, collapse = ",")),
                  r[, "Estimate"])
}
r <- tr("predict upper = 2", predict(fw, newdata = data.frame(x = 0,
                                                              upper = 2),
                                     ndraws = 10))
r <- tr("simulate upper = 2", simulate(fw, newdata = data.frame(x = 0,
                                                                upper = 2),
                                       nsim = 3))
r <- tr("in-sample predict", predict(fw, ndraws = 5))
if (ok(r)) cat("wiener in-sample predict ok\n")
r <- tr("residuals in-sample", residuals(fw))
if (ok(r)) cat("wiener in-sample residuals ok\n")
r <- tr("emmeans", summary(emmeans::emmeans(fw, ~ x, at = list(x = 0),
                                            epred = TRUE)))
if (ok(r)) show("emmeans epred (upper held at its mean?)", r$emmean)
cat("done\n")
