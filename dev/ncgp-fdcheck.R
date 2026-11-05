# The guards dev/ncgp-mut.R no longer kills on the farthest-first factor: how large
# is what gp_field_dtheta() and gp_fd_step() guard against now? A central difference of
# the field at fixed z against the analytic column, on the suite's worst designs,
# and the same with the factor over sorted positions. Seeds in place.
#   Rscript dev/ncgp-fdcheck.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
des <- list(
  grid = local({
    set.seed(17); n <- 120L; xg <- round(runif(n, 0, 10), 1); pos <- sort(unique(xg))
    u <- drop(crossprod(chol(2.25 * exp(-outer(pos, pos, "-")^2 / 8) + diag(1e-8, length(pos))),
                        rnorm(length(pos))))
    data.frame(y = 2 + u[match(xg, pos)] + rnorm(n, 0, 0.5), x = xg)
  }),
  low60 = local({set.seed(20261102); n <- 60L; x <- round(runif(n, 0, 6), 3)
    data.frame(x = x, y = sin(2 * x) + 1e-3 * rnorm(n))}),
  low200 = local({set.seed(3); n <- 200L; x <- runif(n, 0, 6)
    data.frame(x = x, y = sin(2 * x) + 1e-3 * rnorm(n))}),
  long = local({set.seed(4); n <- 300L; x <- runif(n, 0, 10)
    data.frame(x = x, y = 0.3 * x^2 + 0.3 * rnorm(n))})
)
fnm <- getFromNamespace("gp_field", "frmtmb")
dt <- getFromNamespace("gp_field_dtheta", "frmtmb")
ff <- getFromNamespace("gp_factor", "frmtmb")
for (dn in names(des)) {
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = des[[dn]], family = gaussian()))
  bk <- fit$frame$re_blocks[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  z <- fit$estimates$b[bk$b_idx]
  for (ord in c("farthest", "sorted")) {
    b2 <- bk
    if (ord == "sorted") b2$gp_order <- seq_along(z)
    ana <- exp(th[1L]) * ff(th, b2) %*% dt(z, th, b2)
    for (k in seq_along(th)) {
      e <- replace(numeric(length(th)), k, 1e-5)
      fd <- (fnm(z, th + e, b2) - fnm(z, th - e, b2)) / 2e-5
      cat(sprintf("NCGPFD %-7s %-8s theta[%d] max|fd - analytic| %.2e of %.3g\n", dn, ord, k,
                  max(abs(fd - ana[, k])), max(abs(ana[, k]))))
    }
  }
}
