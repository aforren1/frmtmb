# Lane bigx: the optimizer's trace on one fit that ends at "false
# convergence (8)" in the design's coordinates, and the objective's
# noise along the slope's coordinate there.
#   BIGX_ARM=after Rscript dev/bigx-trace2.R <family> <re> <m> <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
a <- commandArgs(TRUE)
fam <- a[1]
re <- as.logical(a[2])
m <- as.numeric(a[3])
s <- as.integer(a[4])
d <- bigx_data(fam, m, s, re = re)
rhs <- if (re) " + (1 | g)" else ""
f <- suppressWarnings(frm(bf(as.formula(paste0("y ~ w", rhs))),
                          family = bigx_family(fam), data = d,
                          control = frmtmb_control(verbose = 2)))
cat("code", f$opt$convergence, f$opt$message, "\n")
B <- f$obj$env$.frmtmb_opt_basis
p <- f$opt$par
q <- drop(B$Mi %*% p)
cat("q at stop:", format(q, digits = 12), "\n")
gq <- drop(crossprod(B$M, as.numeric(f$obj$gr(p))))
cat("gradient in q:", format(gq), "\n")
# the objective along each q coordinate at steps 1e-9 .. 1e-5
for (j in seq_along(q)) {
  for (h in 10^(-9:-5)) {
    qp <- q
    qp[j] <- qp[j] + h
    qm <- q
    qm[j] <- qm[j] - h
    fp <- f$obj$fn(drop(B$M %*% qp))
    fm <- f$obj$fn(drop(B$M %*% qm))
    f0 <- f$obj$fn(drop(B$M %*% q))
    cat(sprintf("  q%d h %.0e  f(+)-f0 %+.3e  f(-)-f0 %+.3e  2nd diff %+.3e\n",
                j, h, fp - f0, fm - f0, (fp + fm - 2 * f0) / h^2))
  }
}
