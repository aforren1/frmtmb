# Reviewer: frm_curve_deriv() reads the extra covariance through the
# stencil, L E L', with E = extra_cov - extra_white formed on the stacked
# rows. The stencil step was tuned for differencing DESIGN rows (1e-6 of
# the range at order 1, 1e-4 at order 2); differencing a covariance
# divides its rounding by e^2 (order 1) or e^4 (order 2). Measured here
# against closed forms on (1) an unseen level of (1 + t | g), whose slope
# variance is Sg[2, 2] and whose second derivative carries none, and (2)
# an exact gp() past its data, whose derivative variances have the
# squared-exponential kernel's closed-form derivatives.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
set.seed(23)
ng <- 25
d <- data.frame(g = factor(rep(seq_len(ng), each = 12)),
                t = rep(seq(0, 1, length.out = 12), ng))
u0 <- rnorm(ng, 0, 0.5); u1 <- rnorm(ng, 0, 0.3)
d$y <- -1 + 2 * d$t + u0[d$g] + u1[d$g] * d$t + rnorm(nrow(d), 0, 0.2)
fit <- frm(bf(y ~ t + (1 + t | g)), data = d)
Sg <- frmtmb::VarCorr(fit)$g$cov[, "Estimate", ]
tt <- seq(0, 1, length.out = 21)
nd1 <- data.frame(t = tt, g = factor("n1", levels = c(levels(d$g), "n1")))
d1 <- frm_curve_deriv(fit, var = "t", newdata = nd1, re_formula = NULL,
                      allow_new_levels = TRUE, simultaneous = FALSE)
d2 <- frm_curve_deriv(fit, var = "t", order = 2, newdata = nd1,
                      re_formula = NULL, allow_new_levels = TRUE,
                      simultaneous = FALSE)
d2s <- frm_curve_deriv(fit, var = "t", order = 2,
                       newdata = data.frame(t = tt, g = factor("3",
                                            levels = levels(d$g))),
                       re_formula = NULL, simultaneous = FALSE)
cat(sprintf(paste0("RE order 1: se/closed form range [%.6f, %.6f] | ",
                   "order 2 at the unseen level: se range [%.3e, %.3e] ",
                   "(closed form 0; the eta's second derivative is %.3e) | ",
                   "order 2 at a seen level: se max %.3e | eps %.3e\n"),
            min(d1$.se / sqrt(vcov(fit)["t", "t"] + Sg[2, 2])),
            max(d1$.se / sqrt(vcov(fit)["t", "t"] + Sg[2, 2])),
            min(d2$.se), max(d2$.se), max(abs(d2$.estimate)),
            max(d2s$.se), attr(d2, "eps")))
# (2) exact gp() past its data
set.seed(5)
x <- sort(runif(50, 0, 4))
dg <- data.frame(x = x, y = sin(1.3 * x) + rnorm(50, 0, 0.15))
fg <- frm(bf(y ~ gp(x)), data = dg)
th <- fg$estimates$theta
s2 <- exp(2 * th[1]); l <- exp(th[2])
pos <- fg$frame$linpreds[["y.mu"]]$gps[[1]]$positions[, 1]
K <- exp(-outer(pos, pos, "-")^2 / (2 * l^2)) + diag(1e-6, length(pos))
xs <- c(3.0, 4.2, 4.8, 5.5)
r <- outer(xs, pos, "-")
k0 <- exp(-r^2 / (2 * l^2))
k1 <- -r / l^2 * k0
k2 <- (r^2 / l^4 - 1 / l^2) * k0
Ki <- solve(K)
# the conditional variance of f' and f'' given the field at the
# positions, the extra part a derivative carries
v1 <- s2 * (1 / l^2 - rowSums((k1 %*% Ki) * k1))
v2 <- s2 * (3 / l^4 - rowSums((k2 %*% Ki) * k2))
lb <- frm_lp_basis(fg, newdata = data.frame(x = pos[1:3]))
for (ord in 1:2) {
  e <- 6 * if (ord == 1) 1e-6 else 1e-4
  for (scale_e in c(1, 10, 100)) {
    ee <- e * scale_e
    st <- data.frame(x = c(xs - ee, xs, xs + ee))
    lbs <- frm_lp_basis(fg, newdata = st, extra_cov = TRUE)
    Es <- lbs$extra_cov - lbs$extra_white
    m <- length(xs)
    L <- matrix(0, m, 3 * m)
    lo <- seq_len(m); mid <- m + lo; hi <- 2 * m + lo
    if (ord == 1) {
      L[cbind(lo, hi)] <- 1 / (2 * ee); L[cbind(lo, lo)] <- -1 / (2 * ee)
    } else {
      L[cbind(lo, hi)] <- 1 / ee^2; L[cbind(lo, mid)] <- -2 / ee^2
      L[cbind(lo, lo)] <- 1 / ee^2
    }
    got <- diag(L %*% Es %*% t(L))
    want <- if (ord == 1) v1 else v2
    cat(sprintf("GP order %d step %.1e: L E L' / closed form %s\n", ord, ee,
                paste(sprintf("%.6g", got / want), collapse = " ")))
  }
}
gd2 <- frm_curve_deriv(fg, var = "x", order = 2,
                       newdata = data.frame(x = xs), simultaneous = FALSE)
gd1 <- frm_curve_deriv(fg, var = "x", order = 1,
                       newdata = data.frame(x = xs), simultaneous = FALSE)
cat("GP frm_curve_deriv order 1 se:", sprintf("%.6g", gd1$.se),
    "| extra part alone sqrt(v1):", sprintf("%.6g", sqrt(v1)), "\n")
cat("GP frm_curve_deriv order 2 se:", sprintf("%.6g", gd2$.se),
    "| extra part alone sqrt(v2):", sprintf("%.6g", sqrt(v2)), "\n")
cat("DONE\n")
