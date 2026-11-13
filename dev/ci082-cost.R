# Lane ci082, fix round: what the fixed-row probe (varcorr_fixed_rows())
# costs in varcorr_trans_rows() on a us() block of dimension 6 (21 theta).
# Two arms in one process, interleaved: the lane's function, and the same
# function with the probe replaced by "nothing fixed". Each round times a
# block of calls long enough to clear the 10 ms clock, 15 rounds, min.
# A control arm (the lane's function against itself) must read 1.0.
#   Rscript dev/ci082-cost.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
set.seed(6)
G <- 60L
d6 <- 6L
S <- 0.4^abs(outer(1:d6, 1:d6, "-"))
U <- matrix(rnorm(G * d6), G) %*% chol(S)
X <- matrix(rnorm(G * 10L * (d6 - 1L)), ncol = d6 - 1L)
colnames(X) <- paste0("x", 1:(d6 - 1L))
dd <- data.frame(g = factor(rep(seq_len(G), each = 10L)), X)
Z <- cbind(1, X)
dd$y <- 1 + rowSums(Z * U[dd$g, ]) + rnorm(nrow(dd), 0, 0.5)
f <- as.formula(paste("y ~ x1 + x2 + x3 + x4 + x5 +",
                      "(1 + x1 + x2 + x3 + x4 + x5 | g)"))
fit <- suppressWarnings(frm(f, family = gaussian(), data = dd))
bk <- fit$frame[["re_blocks"]][[1]]
cat("block", bk[["covstruct"]], "dim", bk[["dim"]], "theta",
    length(bk[["theta_idx"]]), "code", fit$opt$convergence, "\n")
# count the structure's vcov calls in one varcorr_trans_rows()
e <- ns$cs_entry(bk[["covstruct"]])
n_vc <- 0L
vc0 <- e$vcov
cnt <- function(f) function(...) {
  n_vc <<- n_vc + 1L
  f(...)
}
reg <- ns$frmtmb_covstruct_registry$reg
reg[[bk[["covstruct"]]]]$vcov <- cnt(vc0)
assign("reg", reg, envir = ns$frmtmb_covstruct_registry)
invisible(ns$varcorr_trans_rows(fit))
n_with <- n_vc
n_vc <- 0L
invisible(ns$varcorr_fixed_rows(bk, fit$estimates[["theta"]][bk[["theta_idx"]]]))
n_probe <- n_vc
reg[[bk[["covstruct"]]]]$vcov <- vc0
assign("reg", reg, envir = ns$frmtmb_covstruct_registry)
cat("vcov calls per varcorr_trans_rows():", n_with, " of them the probe's:",
    n_probe, "\n")

with_probe <- ns$varcorr_trans_rows
without_probe <- with_probe
environment(without_probe) <- list2env(
  list(varcorr_fixed_rows = function(bk, t0) NULL), parent = ns)
tr_a <- with_probe(fit)
tr_b <- without_probe(fit)
cat("same rows with and without the probe:",
    identical(tr_a, tr_b), "\n")
blk <- function(fn, n) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(n)) fn(fit)
  proc.time()[["elapsed"]] - t0
}
n <- 1L
while (blk(with_probe, n) < 1.2) n <- n * 2L
cat("calls per block:", n, "\n")
res <- matrix(NA_real_, 15, 3, dimnames = list(NULL,
  c("with", "without", "control")))
for (r in 1:15) {
  res[r, "with"] <- blk(with_probe, n)
  res[r, "without"] <- blk(without_probe, n)
  res[r, "control"] <- blk(with_probe, n)
}
m <- apply(res, 2, min)
cat(sprintf("min over 15 rounds, s per call: with %.3e  without %.3e  control %.3e\n",
            m[["with"]] / n, m[["without"]] / n, m[["control"]] / n))
cat(sprintf("ratio with/without %.3f  control/with %.3f\n",
            m[["with"]] / m[["without"]], m[["control"]] / m[["with"]]))
t0 <- proc.time()[["elapsed"]]
for (i in 1:3) invisible(confint_varcorr(fit))
cat(sprintf("confint_varcorr() of the fit: %.3f s per call\n",
            (proc.time()[["elapsed"]] - t0) / 3))
