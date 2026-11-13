# Short check of lane ci082's fix round: the cost of varcorr_fixed_rows()
# inside varcorr_trans_rows() on a us() block of dimension d (6 or 10).
# After dev/ci082-cost.R: interleaved arms in one process, the lane's
# function and the same with the probe replaced by "nothing fixed", plus
# a control (the lane's against itself); blocks of calls of at least
# 1.2 s, 15 rounds, min per arm.
#   Rscript dev/reviews/ci082-short-cost.R <d>
d6 <- as.integer(commandArgs(TRUE)[1])
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
set.seed(6)
G <- if (d6 > 6L) 150L else 60L
per <- if (d6 > 6L) 20L else 10L
S <- 0.4^abs(outer(1:d6, 1:d6, "-"))
U <- matrix(rnorm(G * d6), G) %*% chol(S)
X <- matrix(rnorm(G * per * (d6 - 1L)), ncol = d6 - 1L)
colnames(X) <- paste0("x", 1:(d6 - 1L))
dd <- data.frame(g = factor(rep(seq_len(G), each = per)), X)
Z <- cbind(1, X)
dd$y <- 1 + rowSums(Z * U[dd$g, ]) + rnorm(nrow(dd), 0, 0.5)
xs <- paste(colnames(X), collapse = " + ")
f <- as.formula(paste("y ~", xs, "+ (1 +", xs, "| g)"))
tf <- system.time(fit <- suppressWarnings(frm(f, family = gaussian(),
                                              data = dd)))[["elapsed"]]
bk <- fit$frame[["re_blocks"]][[1]]
cat("block", bk[["covstruct"]], "dim", bk[["dim"]], "theta",
    length(bk[["theta_idx"]]), "code", fit$opt$convergence,
    "fit s", tf, "\n")
t0 <- fit$estimates[["theta"]][bk[["theta_idx"]]]
cat("fixed rows:", sum(ns$varcorr_fixed_rows(bk, t0)), "of",
    length(ns$varcorr_fixed_rows(bk, t0)), "\n")
with_probe <- ns$varcorr_trans_rows
without_probe <- with_probe
environment(without_probe) <- list2env(
  list(varcorr_fixed_rows = function(bk, t0) NULL), parent = ns)
cat("same rows with and without the probe:",
    identical(with_probe(fit), without_probe(fit)), "\n")
blk <- function(fn, n) {
  s <- proc.time()[["elapsed"]]
  for (i in seq_len(n)) fn(fit)
  proc.time()[["elapsed"]] - s
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
cat(sprintf("d %d min over 15 rounds, ms per call: with %.3f without %.3f control %.3f\n",
            d6, 1e3 * m[["with"]] / n, 1e3 * m[["without"]] / n,
            1e3 * m[["control"]] / n))
cat(sprintf("d %d ratio with/without %.3f control/with %.3f\n",
            d6, m[["with"]] / m[["without"]], m[["control"]] / m[["with"]]))
s <- proc.time()[["elapsed"]]
for (i in 1:3) invisible(confint_varcorr(fit))
cat(sprintf("d %d confint_varcorr(): %.4f s per call\n", d6,
            (proc.time()[["elapsed"]] - s) / 3))
