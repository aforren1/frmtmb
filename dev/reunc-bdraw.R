# Lane wt-reunc: the group-effect draw measured DIRECTLY, against
# Q_rr^-1, rather than through the simulated response. A shortfall seen
# in the response variance (dev/reunc-log/plugvar.txt) is either here or
# in the noise draw, and this instrument separates them.
# The library is an ARGUMENT with a default, not a constant: pinned to
# one session's private library this script silently ran against a
# stale build after a later session rebuilt the package elsewhere.
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6"))) 
suppressMessages(library(frmtmb))
set.seed(4)
G <- 6; m <- 4
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 1)[d$g] + rnorm(G * m, 0, 0.8)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
Q <- frmtmb:::joint_precision(fit)
rn <- rownames(Q)
bi <- which(rn == "b")
V <- as.matrix(solve(Q))[bi, bi]
N <- 50000
for (sd in 1:3) {
  set.seed(sd)
  drawer <- frmtmb:::predict_b_drawer(fit, NULL, NULL, N, TRUE)
  B <- t(vapply(seq_len(N), function(s) drawer(fit, s)$estimates[["b"]],
                numeric(length(fit$estimates[["b"]]))))
  Cb <- cov(B)
  # conditional (plug-in) law: Q_bb^-1, not the marginal V above
  Qbb <- as.matrix(Q[bi, bi])
  want <- solve(Qbb)
  cat(sprintf("seed %d plug-in: diag ratio %s\n", sd,
              paste(sprintf("%.4f", diag(Cb) / diag(want)), collapse = " ")))
  cat(sprintf("  mean of draws minus mode (in sd units): %s\n",
              paste(sprintf("%.3f",
                            (colMeans(B) - fit$estimates[["b"]]) /
                              (sqrt(diag(want)) / sqrt(N))),
                    collapse = " ")))
}
