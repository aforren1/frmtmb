# Delayed entry under cox(): the hand likelihood on each knot choice
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.survival)})
set.seed(11)
n <- 200
x <- rnorm(n)
tt <- rexp(n, exp(-0.5 + 0.7 * x))
ct <- rexp(n, 0.2)
d <- data.frame(time = pmin(tt, ct) + 2, cens = as.numeric(tt > ct), x = x)
d$entry <- d$time * 0.5
fit <- frm(bf(time | cens(cens) + trunc(lb = entry) ~ x), family = cox(),
           data = d)
sp_lane <- fit$frame$aterm_values$time$bhaz_spec
r <- range(d$time)
hand <- function(bnd, internal, degree, intercept) {
  I <- function(t) splines2::iSpline(pmin(pmax(t, bnd[1]), bnd[2]),
                                     knots = internal, Boundary.knots = bnd,
                                     degree = degree, intercept = intercept)
  Iy <- I(d$time); Ie <- I(d$entry)
  My <- splines2::mSpline(d$time, knots = internal, Boundary.knots = bnd,
                          degree = degree, intercept = intercept)
  K <- ncol(Iy)
  nll <- function(p) {
    s <- exp(c(0, p[-(1:2)])); s <- s / sum(s)
    eta <- p[1] + p[2] * d$x
    Hy <- exp(eta) * drop(Iy %*% s); He <- exp(eta) * drop(Ie %*% s)
    hy <- exp(eta) * drop(My %*% s)
    -sum(ifelse(d$cens == 0, log(hy), 0) - Hy + He)
  }
  best <- NULL
  for (st in 1:8) {
    set.seed(st)
    o <- optim(c(0, 0, rnorm(K - 1, 0, 0.5)), nll, method = "BFGS",
               control = list(maxit = 10000, reltol = 1e-14))
    if (is.null(best) || o$value < best$value) best <- o
  }
  -best$value
}
cat("lane default knots", sp_lane$boundary, " internal", sp_lane$internal,
    " frmtmb logLik", format(as.numeric(logLik(fit)), digits = 10),
    " hand", format(hand(sp_lane$boundary, sp_lane$internal, sp_lane$degree,
                         sp_lane$intercept), digits = 10), "\n")
bk <- c(max(r[1] - diff(r) / 50, 0), r[2] + diff(r) / 50)
cat("brms default knots", bk, " hand", format(hand(bk, sp_lane$internal,
                                                   sp_lane$degree,
                                                   sp_lane$intercept),
                                              digits = 10), "\n")
