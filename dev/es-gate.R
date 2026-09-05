# P1: the objective gate and expand_b() must agree about a frame that
# predates `has_expand`. Reproduces the review's off-the-mode numbers.
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(frmtmb)
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c)); n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
  }
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
set.seed(42); W <- lattice_W(4, 4); n <- nrow(W)
K <- diag(rowSums(W)) - W + matrix(1 / (1e-3 * n)^2, n, n)
phi <- 1.2 * drop(crossprod(chol(solve(K)), rnorm(n)))
loc <- factor(rep(rownames(W), each = 6), levels = rownames(W))
d <- data.frame(loc = loc, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.5)
fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
           data = d, data2 = list(W = W))

pl <- fit$obj$env$parList(fit$obj$env$last.par.best)
fr <- fit$frame
cat("frame has_expand:", fr[["has_expand"]], " has_rr:", fr[["has_rr"]], "\n")
f_new <- frmtmb:::build_objective(fr)
fr_old <- fr
fr_old[["has_expand"]] <- NULL
cat("stripped has_expand is NULL:", is.null(fr_old[["has_expand"]]), "\n")
f_old <- frmtmb:::build_objective(fr_old)

cat(sprintf("at the mode      : intact %.10f  stripped %.10f  diff %.3e\n",
            f_new(pl), f_old(pl), f_new(pl) - f_old(pl)))
pl2 <- pl
pl2$b <- pl$b + 0.37
cat(sprintf("block shifted .37: intact %.7f  stripped %.7f  diff %.3e\n",
            f_new(pl2), f_old(pl2), f_new(pl2) - f_old(pl2)))
# both flags forced FALSE: frame_needs_expand() still derives TRUE
# from the blocks, which is the robustness the fix buys
bad <- local({
  frx <- fr
  frx[["has_expand"]] <- FALSE
  frx[["has_rr"]] <- FALSE
  frmtmb:::build_objective(frx)
})
cat(sprintf("both flags FALSE, shifted   : %.7f\n", bad(pl2)))
