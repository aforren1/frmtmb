# Case 1 of dev/rel072-v29-diag.R: the theta block of cov.fixed and of
# cov_fixed_prop, se_null, and the jacobian confint_varcorr() uses.
#   Rscript dev/rel072-v29-diag2.R <lib> [<lib2>]
a <- commandArgs(TRUE)
.libPaths(c(a, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(61)
ng <- 12; u <- stats::rnorm(ng, 0, 0.8); per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] + stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
fit <- suppressWarnings(frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
                                 bf(y2 ~ x + (1 | q | id)) + gaussian()),
                            data = dd))
sdr <- frmtmb:::sdr_of(fit)
th <- which(rownames(sdr$cov.fixed) == "theta")
cat("cov.fixed theta block:\n"); print(sdr$cov.fixed[th, th])
cat("cov_fixed_raw theta block:\n"); print(sdr$cov_fixed_raw[th, th])
cat("cov_fixed_prop is null:", is.null(sdr$cov_fixed_prop), "| se_null is null:",
    is.null(sdr$se_null), "\n")
if (!is.null(sdr$cov_fixed_prop)) print(sdr$cov_fixed_prop[th, th])
print(sdr$se_lost)
tr <- frmtmb:::varcorr_trans_rows(fit)
print(tr)
# the jacobian column of the correlation's parameter, as
# varcorr_trans_rows() forms it, times the step: what each row moved
bk <- fit$frame$re_blocks[[1]]
t0 <- fit$estimates$theta[bk$theta_idx]
gfun <- function(tt) {
  V <- frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(tt, bk)
  C <- stats::cov2cor(V)
  c(log(sqrt(diag(V))), atanh(pmin(pmax(C[lower.tri(C)], -0.9999), 0.9999)))
}
h <- 1e-5 * max(abs(t0[3]), 1)
tp <- t0; tp[3] <- tp[3] + h; tm <- t0; tm[3] <- tm[3] - h
cat("moved along theta_3 over the step 2h:", format(gfun(tp) - gfun(tm), digits = 3),
    "| 16 ulps of g0:", format(16 * .Machine$double.eps * pmax(abs(gfun(t0)), 1), digits = 3), "\n")
