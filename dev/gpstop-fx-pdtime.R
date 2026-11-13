# Lane gpstop, short check: the cost of the inner-Hessian test at q random
# effects, the sparse Cholesky of H - 1e-2 I (gp_inner_pd(), shipped)
# against the dense eigendecomposition it replaced, in one process, arms
# interleaved, 3 rounds, minimum taken. Model: poisson y ~ gp(x) + (1 | s)
# with 60 positions and q - 60 subjects, one row each, built at its start
# (iter.max = 0); set.seed(9400 + q). Args: q values (comma).
.libPaths(c("C:/Users/adf44/source/r/wt-gpstop-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
qs <- as.integer(strsplit(commandArgs(TRUE)[1], ",")[[1]])
ns <- asNamespace("frmtmb")
dense <- function(obj, par) {
  env <- obj$env
  obj$fn(par)
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  min(eigen(0.5 * (H + t(H)), symmetric = TRUE, only.values = TRUE)$values)
}
for (q in qs) {
  set.seed(9400 + q)
  m <- q - 60L
  xs <- seq(0, 10, length.out = 60)
  d <- data.frame(x = xs[sample(60, m, TRUE)], s = factor(seq_len(m)))
  d$y <- rpois(m, exp(1 + sin(d$x) + rnorm(m, sd = 0.3)))
  fit <- suppressWarnings(frm(y ~ gp(x) + (1 | s), data = d,
                              family = poisson(),
                              control = frmtmb_control(
                                optCtrl = list(iter.max = 0))))
  obj <- fit$obj
  par <- fit$opt$par
  nr <- length(obj$env$random)
  best <- c(sparse = Inf, dense = Inf)
  mem <- c(sparse = NA, dense = NA)
  out <- list()
  for (r in 1:3) for (arm in if (r %% 2) c("sparse", "dense") else
                                c("dense", "sparse")) {
    invisible(gc(reset = TRUE))
    t0 <- proc.time()[[3]]
    v <- if (arm == "sparse") ns$gp_inner_pd(obj, par, 1e-2) else
      dense(obj, par)
    el <- proc.time()[[3]] - t0
    g <- gc()
    best[arm] <- min(best[arm], el)
    mem[arm] <- sum(g[, ncol(g)])
    out[[arm]] <- v
  }
  cat(sprintf(paste0("PDTIME q %d random %d | sparse %.3f s (max MB %.0f) ",
                     "verdict %s | dense %.3f s (max MB %.0f) ev_min %.4g ",
                     "agree %s\n"),
              q, nr, best[["sparse"]], mem[["sparse"]], out$sparse,
              best[["dense"]], mem[["dense"]], out$dense,
              identical(out$sparse, out$dense > 1e-2)))
  # both sides of the answer: cuts just below and just above the eigenvalue
  cat(sprintf("PDEDGE q %d cut 0.99 ev %s (want TRUE) | cut 1.01 ev %s (want FALSE)\n",
              q, ns$gp_inner_pd(obj, par, 0.99 * out$dense),
              ns$gp_inner_pd(obj, par, 1.01 * out$dense)))
  flush.console()
}
cat("DONE\n")
