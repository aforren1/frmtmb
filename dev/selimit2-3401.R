# Lane selimit2, backlog item 3401 (rerun on 0.82.0): lane adhess's nl_re
# model (dev/adhess-restep.R, seed 11), a nonlinear body with a random
# effect. The fit's own outer standard errors against a Richardson
# Jacobian of the AD gradient (numDeriv, six halvings, base step 1e-2 in
# standard-error units), as dev/adhess-restep.R's reference.
#   Rscript dev/selimit2-3401.R <lib or base>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(11)
n <- 300
d <- data.frame(x = 5 + runif(n), g = factor(rep(1:15, length.out = n)))
d$y <- (2 + rnorm(15, 0, 0.2)[d$g]) * exp(0.8 * (d$x - 5)) +
  rnorm(n, 0, 0.3)
f <- suppressWarnings(frm(bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1,
                             nl = TRUE), family = gaussian(), data = d,
                          start = list(beta = c(2 * exp(-0.8 * 5), 0.8))))
se_fit <- sqrt(diag(vcov(f, full = TRUE)))
obj <- f$obj
p <- ns$sdr_outer_point(obj)
saved <- ns$obj_state_save(obj)
H0 <- ns$fd_hessian(p, function(q) obj$gr(q))$H
D <- sqrt(abs(diag(H0)))
ref <- function(d0) {
  g0 <- function(dz) as.numeric(obj$gr(p + dz / D)) / D
  J <- numDeriv::jacobian(g0, rep(0, length(p)), method = "Richardson",
                          method.args = list(eps = d0, d = d0, zero.tol = 1,
                                             r = 6, v = 2))
  J <- (J + t(J)) / 2
  J * outer(D, D)
}
R1 <- ref(1e-1)
R2 <- ref(1e-2)
ns$obj_state_restore(obj, saved)
s1 <- sqrt(diag(solve(R1)))
s2 <- sqrt(diag(solve(R2)))
cat("par     ", names(p), "\n")
cat("se fit  ", sprintf("%.9g", se_fit), "\n")
cat("se ref  ", sprintf("%.9g", s2), "\n")
cat("refself max rel", max(abs(s1 / s2 - 1)), "\n")
cat("fit vs ref max rel", max(abs(se_fit / s2 - 1)), "\n")
