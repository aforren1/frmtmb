# Lane bigx: what opt_basis() costs a fit whose design is well
# conditioned (it runs reparam_basis()'s QR and returns NULL), measured
# on its own and as a share of the fit. Interleaved rounds, the minimum
# of each; a control (the same call twice) must read 1.
#   BIGX_ARM=after Rscript dev/bigx-cost.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
set.seed(1)
for (shape in list(c(n = 1e4, p = 20), c(n = 1e5, p = 50))) {
  n <- shape[["n"]]
  p <- shape[["p"]]
  X <- matrix(rnorm(n * (p - 1)), n, p - 1)
  colnames(X) <- paste0("x", seq_len(p - 1))
  d <- data.frame(X, y = rpois(n, exp(0.2 + X[, 1] * 0.1)))
  fo <- as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
  f <- frm(bf(fo), family = poisson(), data = d)
  tm <- function(expr, k) {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(k)) force(expr())
    (proc.time()[["elapsed"]] - t0) / k
  }
  k <- if (n > 5e4) 3L else 20L
  basis <- function() ns$opt_basis(f$frame, f$obj, "nlminb")
  rounds <- replicate(5, c(basis = tm(basis, k), control = tm(basis, k)))
  fit_t <- min(replicate(3, system.time(
    frm(bf(fo), family = poisson(), data = d))[["elapsed"]]))
  cat(sprintf(paste0("n %g p %g: opt_basis() %.4f s (control %.4f s, ",
                     "ratio %.3f), basis NULL %s; whole fit %.3f s; ",
                     "share %.1f%%\n"),
              n, p, min(rounds["basis", ]), min(rounds["control", ]),
              min(rounds["basis", ]) / min(rounds["control", ]),
              is.null(basis()), fit_t, 100 * min(rounds["basis", ]) / fit_t))
}
