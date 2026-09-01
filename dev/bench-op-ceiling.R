## Theoretical ceiling for overlapping fn and gr.
##
## TMB's obj$gr(x) needs the inner Newton solve at x. Sequentially,
## obj$fn(x) has just done it and gr reuses it (last.par cache), so the
## pair costs solve + adjoint. In optimParallel the two land on separate
## processes with separate caches, so the gr worker pays solve + adjoint
## on its own and the round costs max(solve, solve + adjoint). Measure
## the three quantities and the implied ceiling.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
suppressPackageStartupMessages(library(frmtmb))
d <- lme4::InstEval
form <- y ~ service + (1 | s) + (1 | d)

fit <- suppressWarnings(frm(form, data = d, family = gaussian(),
                            se = FALSE,
                            control = frmtmb_control(restarts = 0)))
obj <- fit$obj
p <- fit$opt$par
n <- 20L
set.seed(11)
jit <- lapply(seq_len(n), function(i) p + rnorm(length(p), 0, 0.02))

tm <- function(expr) {
  t <- proc.time()[["elapsed"]]; force(expr)
  proc.time()[["elapsed"]] - t
}

t_fn <- t_gr_cached <- t_gr_fresh <- numeric(n)
for (i in seq_len(n)) {
  q <- jit[[i]]
  t_fn[i] <- tm(obj$fn(q))            # inner Newton solve at a new point
  t_gr_cached[i] <- tm(obj$gr(q))     # adjoint only (last.par == q)
}
for (i in seq_len(n)) {               # fresh point, gr with no prior fn
  q <- jit[[i]] + rnorm(length(p), 0, 0.02)
  t_gr_fresh[i] <- tm(obj$gr(q))
}
m_fn <- median(t_fn); m_grc <- median(t_gr_cached)
m_grf <- median(t_gr_fresh)
cat(sprintf("fn(x)  [inner solve]            : %.4f s\n", m_fn))
cat(sprintf("gr(x)  after fn(x) [adjoint]    : %.4f s\n", m_grc))
cat(sprintf("gr(x)  at a fresh x [solve+adj] : %.4f s\n", m_grf))
cat(sprintf("\nsequential pair (fn then gr)   : %.4f s\n", m_fn + m_grc))
cat(sprintf("parallel round max(fn, gr_fresh): %.4f s\n",
            max(m_fn, m_grf)))
cat(sprintf("ceiling speedup on evaluations : %.2fx\n",
            (m_fn + m_grc) / max(m_fn, m_grf)))
cat(sprintf("\n(if gr were free of the solve, the ceiling would be %.2fx)\n",
            (m_fn + m_grc) / max(m_fn, m_grc)))

## round-trip cost of the parallel dispatch itself: what optimParallel
## pays per evalFG on top of the evaluations
suppressPackageStartupMessages(library(parallel))
cl <- makePSOCKcluster(2)
invisible(clusterEvalQ(cl, NULL))
rt <- replicate(50, tm(clusterEvalQ(cl, 1)))
cat(sprintf("\nempty 2-worker clusterEvalQ round trip: median %.4f s\n",
            median(rt)))
x <- rnorm(5)
rt2 <- replicate(50, tm({
  clusterExport(cl, "x", environment())
  clusterEvalQ(cl, 1)
}))
cat(sprintf("export a 5-vector + round trip        : median %.4f s\n",
            median(rt2)))
stopCluster(cl)
