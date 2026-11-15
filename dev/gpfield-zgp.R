# Lane gpfield, item 5319: an exact gp()'s z in frmtmb against brms 2.23.0's
# zgp for the same field. brms's z is L_b^-1 f / sdgp with L_b the Cholesky
# factor of the kernel over brms's positions in brms's order (standata's
# Xgp_1, first appearance); frmtmb's is fit$estimates$b, z over its
# farthest-first order. Same field f, read from frmtmb's fit. Data as
# test-brms-lik-gp-exact.R's gpx_data() (set.seed(20261102), 80 rows on 13
# integer positions, unsorted). Args: none
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
suppressPackageStartupMessages(library(brms))
set.seed(20261102)
n <- 80
d <- data.frame(x = sample(0:12, n, TRUE), z = sample(0:12, n, TRUE))
d$y <- 2 * sin(0.7 * d$x) + stats::rnorm(n, 0, 0.3)
fit <- frmtmb::frm(frmtmb::bf(y ~ gp(x)), data = d, family = gaussian())
ns <- asNamespace("frmtmb")
bk <- fit$frame$re_blocks[[1]]
th <- fit$estimates$theta[bk$theta_idx]
z_f <- fit$estimates$b[bk$b_idx]
f <- ns$gp_field(z_f, th, bk)        # at frmtmb's sorted positions
pos <- fit$frame$linpreds[[1]]$gps[[1]]$positions[, 1]
sd <- brms::standata(brms::bf(y ~ gp(x)), data = d, family = gaussian())
xb <- as.numeric(sd$Xgp_1) * bk$gp_lscale_div
io <- match(round(xb, 8), round(pos, 8))
K <- ns$gp_corr(th, bk)
Lb <- t(chol(K[io, io]))
z_b <- as.numeric(forwardsolve(Lb, f[io])) / exp(th[1])
cat(sprintf("ZGP positions %d | frmtmb order of the positions %s\n",
            length(pos), paste(bk$gp_order, collapse = ",")))
cat(sprintf("ZGP brms order of the positions (sorted index) %s\n",
            paste(io, collapse = ",")))
for (i in seq_along(z_f)) {
  cat(sprintf("ZGP i %2d | frmtmb zgp[i] %9.5f | brms zgp[i] %9.5f\n", i,
              z_f[i], z_b[i]))
}
f_b <- exp(th[1]) * as.numeric(Lb %*% z_b)
cat(sprintf(paste0("ZGP max |frmtmb zgp - brms zgp| %.4f | the two fields' ",
                   "max difference %.2e | sum of squares %.6f %.6f\n"),
            max(abs(z_f - z_b)), max(abs(f_b - f[io])), sum(z_f^2),
            sum(z_b^2)))
