# Short check, MAJOR 1: delayed entry and trunc(ub) under cox().
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.survival)})
cat("survival from", find.package("frmtmb.survival"), "\n")
set.seed(11)
n <- 200
x <- rnorm(n)
tt <- rexp(n, exp(-0.5 + 0.7 * x))
ct <- rexp(n, 0.2)
d <- data.frame(time = pmin(tt, ct) + 2, cens = as.numeric(tt > ct), x = x)
d$entry <- d$time * 0.5
r <- range(d$time)
bk <- c(max(r[1] - diff(r) / 50, 0), r[2] + diff(r) / 50)
fit <- frm(bf(time | cens(cens) + trunc(lb = entry) ~ x), family = cox(),
           data = d)
sp <- fit$frame$aterm_values$time$bhaz_spec
cat("DEFAULT boundary", format(sp$boundary, digits = 10), " brms",
    format(bk, digits = 10), " identical", identical(sp$boundary, bk),
    " logLik", format(as.numeric(logLik(fit)), digits = 10), "\n")
fg <- tryCatch(frm(bf(time | cens(cens) + trunc(lb = entry) ~ x),
                   family = cox(Boundary.knots = bk), data = d),
               error = function(e) conditionMessage(e))
if (is.character(fg)) cat("GIVEN brms knots ERROR", fg, "\n") else
  cat("GIVEN brms knots accepted, logLik",
      format(as.numeric(logLik(fg)), digits = 10), " identical to default",
      identical(as.numeric(logLik(fg)), as.numeric(logLik(fit))), "\n")
# hand likelihood with weights that CAN reach 0: w >= 0 by L-BFGS-B
bnd <- sp$boundary
I <- function(t) splines2::iSpline(pmin(pmax(t, bnd[1]), bnd[2]),
                                   knots = sp$internal, Boundary.knots = bnd,
                                   degree = sp$degree,
                                   intercept = sp$intercept)
Iy <- I(d$time); Ie <- I(d$entry)
My <- splines2::mSpline(d$time, knots = sp$internal, Boundary.knots = bnd,
                        degree = sp$degree, intercept = sp$intercept)
K <- ncol(Iy)
nll_w <- function(p) {
  w <- p[-(1:2)]
  if (sum(w) <= 0) return(1e10)
  s <- w / sum(w)
  eta <- p[1] + p[2] * d$x
  Hy <- exp(eta) * drop(Iy %*% s); He <- exp(eta) * drop(Ie %*% s)
  hy <- exp(eta) * drop(My %*% s)
  v <- -sum(ifelse(d$cens == 0, log(hy), 0) - Hy + He)
  if (!is.finite(v)) 1e10 else v
}
best <- NULL
for (st in 1:12) {
  set.seed(st)
  o <- optim(c(0, 0, runif(K, 0.1, 1)), nll_w, method = "L-BFGS-B",
             lower = c(-Inf, -Inf, rep(0, K)),
             control = list(maxit = 10000, factr = 1e3))
  if (is.null(best) || o$value < best$value) best <- o
}
s_hat <- best$par[-(1:2)] / sum(best$par[-(1:2)])
cat("HAND L-BFGS-B (weights may be 0) max logLik",
    format(-best$value, digits = 10), " weights", format(s_hat, digits = 3),
    " beta_x", format(best$par[2], digits = 6), " frmtmb beta_x",
    format(fixef(fit)["x", "Estimate"], digits = 6), "\n")
## trunc(ub) past the data: default widens, short given knots refused
dU <- d
dU$ub <- dU$time * 1.6
fu <- frm(bf(time | cens(cens) + trunc(ub = ub) ~ x), family = cox(),
          data = dU)
cat("TRUNC UB default boundary", format(fu$frame$aterm_values$time$
                                          bhaz_spec$boundary, digits = 8),
    " max ub", format(max(dU$ub), digits = 8), " brms upper",
    format(bk[2], digits = 8), "\n")
ru <- tryCatch({frm(bf(time | cens(cens) + trunc(ub = ub) ~ x),
                    family = cox(Boundary.knots = bk), data = dU); "accepted"},
               error = function(e) substr(conditionMessage(e), 1, 110))
cat("TRUNC UB given brms knots:", ru, "\n")
## both: delayed entry and an upper bound
dB <- dU
fb <- frm(bf(time | cens(cens) + trunc(lb = entry, ub = ub) ~ x),
          family = cox(), data = dB)
cat("BOTH default boundary", format(fb$frame$aterm_values$time$bhaz_spec$
                                      boundary, digits = 8), "\n")
