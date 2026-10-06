# Reviewer round 2, continuation of dev/gpby-rev2-deriv.R: difference
# curves of derivatives at one x grid (frmtmb.spline refuses two
# different x grids by design): two unseen levels, and a seen level
# against an unseen one, orders 1 and 2, .se^2 against D V D' plus the
# extra part written by hand.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
set.seed(6)
ng <- 20
d4 <- data.frame(g = factor(rep(seq_len(ng), each = 8)),
                 x = round(runif(8 * ng, 0, 4), 2))
u0 <- rnorm(ng, 0, 0.6); u1 <- rnorm(ng, 0, 0.3)
d4$y <- sin(d4$x) + u0[d4$g] + u1[d4$g] * d4$x + rnorm(nrow(d4), 0, 0.2)
f4 <- frm(bf(y ~ gp(x) + (1 + x | g)), data = d4)
Sg <- frmtmb::VarCorr(f4)$g$cov[, "Estimate", ]
lev <- c(levels(d4$g), "n1", "n2")
xs <- c(1.111, 4.5, 5.3)
mk <- function(l) data.frame(x = xs, g = factor(l, levels = lev))
Dref <- function(a, b, e, o) {
  A <- function(h) as.matrix(frm_lp_basis(f4, newdata = h,
                                          allow_new_levels = TRUE)$A)
  st <- function(g) {
    lo <- g; hi <- g; lo$x <- g$x - e; hi$x <- g$x + e
    if (o == 1) (A(hi) - A(lo)) / (2 * e) else
      (A(hi) - 2 * A(g) + A(lo)) / e^2
  }
  D <- st(a) - st(b)
  diag(D %*% frm_lp_basis(f4, newdata = a, allow_new_levels = TRUE)$V %*%
         t(D))
}
for (o in 1:2) {
  e <- if (o == 1) 4.2e-6 else 4.2e-4
  for (cs in list(c("n1", "n2", 2), c("3", "n1", 1))) {
    dv <- frm_curve_deriv(f4, var = "x", order = o, newdata = mk(cs[1]),
                          contrast = mk(cs[2]), re_formula = NULL,
                          allow_new_levels = TRUE, simultaneous = FALSE,
                          eps = e)
    want <- Dref(mk(cs[1]), mk(cs[2]), e, o) +
      as.numeric(cs[3]) * (if (o == 1) Sg[2, 2] else 0)
    cat(sprintf(paste0("deriv difference %s - %s, order %d: .se^2 / ",
                       "(D V D' + k Sg22) %s | .se %s\n"), cs[1], cs[2], o,
                paste(sprintf("%.8f", dv$.se^2 / want), collapse = " "),
                paste(sprintf("%.3e", dv$.se), collapse = " ")))
  }
}
cat("DONE\n")
