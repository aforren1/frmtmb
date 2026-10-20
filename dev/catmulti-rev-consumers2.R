# Reviewer: the multinomial clamp near p = 0, the multinomial routes the
# first script did not reach, coef() on a categorical fit, and
# frmtmb.sample draws under refcat = "b" and NA. Lane library.
# Log dev/catmulti-rev-out/consumers2.txt.
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb from", dirname(find.package("frmtmb")), "sample from",
    dirname(find.package("frmtmb.sample")), "\n")
source("dev/catmulti-rev-data.R")
sec <- function(t) cat("\n========", t, "========\n")
shape <- function(nm, x) {
  if (is.null(dim(x))) { cat(nm, ": length", length(x), "class", class(x)[1], "\n"); return(invisible()) }
  cat(nm, ": dim", paste(dim(x), collapse = " x "), "\n")
  dn <- dimnames(x)
  for (k in seq_along(dn)) cat("  dimnames[[", k, "]]: ",
    if (is.null(dn[[k]])) "NULL" else paste(utils::head(dn[[k]], 8), collapse = ", "), "\n", sep = "")
}
TRY <- function(nm, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err")),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") },
    message = function(cm) invokeRestart("muffleMessage"))
  if (inherits(r, "err")) cat(nm, ": ERROR:", r, "\n") else if (is.list(r) && !is.data.frame(r))
    cat(nm, ": list", names(r), "\n") else shape(nm, r)
  if (length(w)) cat("  WARNINGS:", substr(unique(w), 1, 200), sep = "\n    ")
  invisible(if (inherits(r, "err")) NULL else r)
}

sec("multinomial, a rare category, small trials")
set.seed(41)
n <- 300
dr <- data.frame(x = rnorm(n))
E <- cbind(0, 0.2 + 0.5 * dr$x, -3.5 - 1.5 * dr$x)
P <- exp(E) / rowSums(exp(E))
dr$N <- sample(1:6, n, TRUE)
Yc <- t(vapply(seq_len(n), function(i) as.numeric(rmultinom(1, dr$N[i], P[i, ])), numeric(3)))
colnames(Yc) <- c("y1", "y2", "y3"); dr$Yc <- Yc
cat("column totals:", colSums(Yc), "\n")
fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dr)
print(round(fixef(fm), 4))
ft <- fitted(fm)
est <- ft[, "Estimate", ]; se <- ft[, "Est.Error", ]
wl <- est - qnorm(0.975) * se; wh <- est + qnorm(0.975) * se
lo <- ft[, "Q2.5", ]; hi <- ft[, "Q97.5", ]
cat("Wald lower < 0 in", sum(wl < 0), "of", length(wl), "cells; Wald upper > trials in",
    sum(wh > dr$N), "\n")
cat("reported lo == 0 in", sum(lo == 0), "; hi == trials in", sum(hi == dr$N), "\n")
cat("reported lo > est - 1.96 se (clamped) cells:", sum(lo > wl + 1e-12), "\n")
k <- which(wl < 0, arr.ind = TRUE)
for (r in head(seq_len(nrow(k)), 4)) {
  i <- k[r, 1]; j <- k[r, 2]
  cat(sprintf("row %d cat %d N %d x %.2f est %.5f se %.5f Wald [%.5f, %.5f] reported [%.5f, %.5f]\n",
              i, j, dr$N[i], dr$x[i], est[i, j], se[i, j], wl[i, j], wh[i, j], lo[i, j], hi[i, j]))
}
# a probability from categorical(): same check, pre-existing clamp
TRY("fitted linear mult", fitted(fm, scale = "linear"))
pr <- TRY("predict mult", predict(fm))
TRY("fitted newdata no trials", fitted(fm, newdata = data.frame(x = 0)))
fn <- TRY("fitted newdata trials 5", fitted(fm, newdata = data.frame(x = c(0, 3), N = c(5, 1))))
if (!is.null(fn)) print(round(fn[, , 3], 5))
TRY("dharma", dharma_residuals(fm))
ce <- TRY("ce mult", conditional_effects(fm, "x", categorical = TRUE)[[1]])
if (!is.null(ce)) { cat("ce cols:", names(ce), "\n"); print(head(ce[, intersect(names(ce),
  c("x", "cats__", "N", "estimate__", "se__", "lower__", "upper__"))], 3)) }
TRY("ce mult method predict", conditional_effects(fm, "x", categorical = TRUE, method = "predict"))
TRY("ce mult no categorical", conditional_effects(fm, "x"))

sec("coef() on categorical with (1 | p | g)")
d3 <- rev_cat_data(21, 3)
fb <- frm(bf(y ~ x + (1 | p | g)), family = categorical(refcat = "b"), data = d3)
cf <- coef(fb)
cat("coef names:", names(cf), "\n")
str(cf, max.level = 2, give.attr = FALSE)
TRY("ranef", ranef(fb))

sec("frmtmb.sample draws, refcat = 'b' and NA, no group term")
dc <- rev_cat_data(23, 3, n = 300)
fb0 <- frm(bf(y ~ x), family = categorical(refcat = "b"), data = dc)
fa0 <- frm(bf(y ~ x), family = categorical(), data = dc)
pri <- c(set_prior("normal(0, 2)", class = "Intercept", dpar = "mua"),
         set_prior("normal(0, 2)", class = "Intercept", dpar = "mub"),
         set_prior("normal(0, 2)", class = "Intercept", dpar = "muc"),
         set_prior("normal(0, 2)", class = "b", dpar = "mua"),
         set_prior("normal(0, 2)", class = "b", dpar = "mub"),
         set_prior("normal(0, 2)", class = "b", dpar = "muc"))
fn0 <- frm(bf(y ~ x), family = categorical(refcat = NA), data = dc, prior = pri)
smp <- function(f) suppressWarnings(suppressMessages(
  frm_sample(f, chains = 1, iter = 400, warmup = 200, refresh = 0, seed = 9)))
for (nm in c("refb", "refa", "refNA")) {
  f <- switch(nm, refb = fb0, refa = fa0, refNA = fn0)
  ds <- TRY(paste(nm, "frm_sample"), smp(f))
  if (is.null(ds)) next
  cat(nm, "variables:", head(grep("^b_", variables(ds), value = TRUE), 8), "\n")
  ep <- TRY(paste(nm, "posterior_epred"), posterior_epred(ds))
  pp <- TRY(paste(nm, "posterior_predict"), posterior_predict(ds))
  if (!is.null(pp)) cat("  posterior_predict values:", head(sort(unique(as.vector(pp)))), "\n")
  ll <- TRY(paste(nm, "log_lik"), log_lik(ds))
  pl <- TRY(paste(nm, "posterior_linpred"), posterior_linpred(ds))
  TRY(paste(nm, "linpred transform"), posterior_linpred(ds, transform = TRUE))
  TRY(paste(nm, "linpred dpar mua"), posterior_linpred(ds, dpar = "mua"))
  if (!is.null(ep)) {
    Pm <- apply(ep, c(2, 3), mean)
    cat("  max |mean epred - fitted|:", max(abs(Pm - fitted(f)[, "Estimate", ])), "\n")
    cat("  epred cats dimnames:", dimnames(ep)[[3]], "\n")
  }
  if (!is.null(ll) && !is.null(ep)) {
    # log_lik must be the log of epred at the observed category, per draw
    yi <- as.integer(dc$y)
    lp <- log(sapply(seq_len(nrow(dc)), function(i) ep[, i, yi[i]]))
    cat("  max |log_lik - log epred[y]|:", max(abs(ll - lp)), "\n")
  }
  if (!is.null(pp)) {
    cat("  posterior_predict shares:", round(prop.table(table(factor(pp, levels = 1:3))), 4),
        " epred mean:", round(colMeans(apply(ep, c(2, 3), mean)), 4), "\n")
  }
}
sec("frmtmb.sample draws: two-category categorical, transform = TRUE")
d2 <- dc; d2$y <- factor(ifelse(d2$y == "a", "a", "b"))
f2 <- frm(bf(y ~ x), family = categorical(), data = d2)
ds2 <- TRY("2cat frm_sample", smp(f2))
if (!is.null(ds2)) {
  TRY("2cat posterior_linpred", posterior_linpred(ds2))
  r <- TRY("2cat linpred transform", posterior_linpred(ds2, transform = TRUE))
  if (!is.null(r)) cat("  range:", range(r), "\n")
}
sec("frmtmb.sample: mixture posterior_linpred, transform")
set.seed(8)
dmx <- data.frame(x = rnorm(300))
dmx$y <- ifelse(runif(300) < 0.4, rnorm(300, -2 + dmx$x), rnorm(300, 2 + 0.5 * dmx$x))
fmx <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = dmx)
dsm <- TRY("mix frm_sample", smp(fmx))
if (!is.null(dsm)) {
  TRY("mix posterior_linpred", posterior_linpred(dsm))
  TRY("mix linpred transform", posterior_linpred(dsm, transform = TRUE))
}
TRY("mix fitted linear", fitted(fmx, scale = "linear"))
cat("DONE\n")
