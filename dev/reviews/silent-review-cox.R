# Review of item 3874. Usage: Rscript silent-review-cox.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.survival)})
cat("arm", arm, "survival from", find.package("frmtmb.survival"), "\n")

bounds_data <- function(seed = 11, n = 200) {
  set.seed(seed)
  x <- stats::rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  d <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                  tt = tt)
  d$hi <- d$time
  top <- order(d$time, decreasing = TRUE)[1:2]
  d$cens[top] <- 2
  d$hi[top] <- max(d$time) * c(1.5, 2)
  mid <- setdiff(seq(5, n, by = 10), top)
  d$cens[mid] <- 2
  d$hi[mid] <- d$time[mid] * 1.2
  d
}

## 1. no extra bound: frmtmb's Zbhaz / Zcbhaz against brms's standata
d0 <- bounds_data()
d0$cens[d0$cens == 2] <- 0
fr0 <- frm(bf(time | cens(cens) ~ x), family = cox(), data = d0,
           dry_run = "frame")
av0 <- fr0$aterm_values$time
sd0 <- brms::make_standata(brms::bf(time | cens(cens) ~ x),
                           family = brms::cox(), data = d0)
cat("brms standata names with bhaz:",
    grep("bhaz", names(sd0), value = TRUE), "\n")
zb <- sd0[[grep("^Zbhaz", names(sd0), value = TRUE)[1]]]
cat("NOBOUND Zbhaz identical:", identical(unname(as.matrix(av0$Zbhaz)),
                                          unname(as.matrix(zb))),
    " max abs diff", max(abs(unname(as.matrix(av0$Zbhaz)) -
                               unname(as.matrix(zb)))), "\n")
cat("NOBOUND boundary", av0$bhaz_spec$boundary, " internal",
    av0$bhaz_spec$internal, "\n")

## 2. interval bound past the data: the fit against a hand likelihood on
## the same knots, and the baseline survivor against the truth
d <- bounds_data()
fit <- frm(bf(time | cens(cens, hi) ~ x), family = cox(), data = d)
sp <- fit$frame$aterm_values$time$bhaz_spec
cat("FIT boundary", sp$boundary, " internal", sp$internal, " logLik",
    format(as.numeric(logLik(fit)), digits = 10), "\n")
I <- function(t) {
  tb <- pmin(pmax(t, sp$boundary[1]), sp$boundary[2])
  splines2::iSpline(tb, knots = sp$internal, Boundary.knots = sp$boundary,
                    degree = sp$degree, intercept = sp$intercept)
}
M <- function(t) splines2::mSpline(t, knots = sp$internal,
                                   Boundary.knots = sp$boundary,
                                   degree = sp$degree,
                                   intercept = sp$intercept)
Iy <- I(d$time); My <- M(d$time); I2 <- I(d$hi)
K <- ncol(Iy)
nll <- function(p) {
  b0 <- p[1]; b1 <- p[2]
  s <- exp(c(0, p[-(1:2)])); s <- s / sum(s)
  eta <- b0 + b1 * d$x
  Hy <- exp(eta) * drop(Iy %*% s)
  H2 <- exp(eta) * drop(I2 %*% s)
  hy <- exp(eta) * drop(My %*% s)
  l <- ifelse(d$cens == 0, log(hy) - Hy,
       ifelse(d$cens == 1, -Hy,
       ifelse(d$cens == -1, log1p(-exp(-Hy)),
              -Hy + log1p(-exp(-(H2 - Hy))))))
  -sum(l)
}
best <- NULL
for (st in 1:5) {
  set.seed(st)
  o <- optim(c(0, 0, rnorm(K - 1, 0, 0.3)), nll, method = "BFGS",
             control = list(maxit = 5000, reltol = 1e-14))
  if (is.null(best) || o$value < best$value) best <- o
}
cat("HAND max logLik", format(-best$value, digits = 10), " beta_x",
    format(best$par[2], digits = 6), " frmtmb beta_x",
    format(fixef(fit)["x", "Estimate"], digits = 6), "\n")
# baseline survivor at the last bounds vs the truth exp(-exp(-0.5) t)
s <- exp(c(0, best$par[-(1:2)])); s <- s / sum(s)
tq <- c(max(d$time), max(d$hi) * 0.75, max(d$hi))
S_hat <- exp(-exp(best$par[1]) * drop(I(tq) %*% s))
cat("S0 at", format(tq, digits = 4), ": fitted", format(S_hat, digits = 4),
    " true", format(exp(-exp(-0.5) * tq), digits = 4), "\n")

## 3. delayed entry below brms's lower knot: the default knots move off
## brms's, and given knots at brms's default are refused
dt <- d0
dt$time <- dt$time + 2
dt$entry <- dt$time * 0.5
fT <- tryCatch(frm(bf(time | cens(cens) + trunc(lb = entry) ~ x),
                   family = cox(), data = dt),
               error = function(e) conditionMessage(e))
if (is.character(fT)) cat("TRUNC default fit ERROR:", fT, "\n") else
  cat("TRUNC default boundary", fT$frame$aterm_values$time$bhaz_spec$boundary,
      " logLik", format(as.numeric(logLik(fT)), digits = 10), "\n")
r <- range(dt$time)
bk <- c(max(r[1] - diff(r) / 50, 0), r[2] + diff(r) / 50)
fB <- tryCatch(frm(bf(time | cens(cens) + trunc(lb = entry) ~ x),
                   family = cox(Boundary.knots = bk), data = dt),
               error = function(e) conditionMessage(e))
if (is.character(fB)) cat("TRUNC brms knots given: ERROR:",
                          substr(fB, 1, 200), "\n") else
  cat("TRUNC brms knots given: boundary",
      fB$frame$aterm_values$time$bhaz_spec$boundary, " logLik",
      format(as.numeric(logLik(fB)), digits = 10), "\n")
