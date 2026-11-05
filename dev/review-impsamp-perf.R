# Item 8: cost of one freeze and one evaluation on large designs, interleaved
# arms in separate processes (dev/review-impsamp-perf.sh), with a fixed
# arithmetic control. G x n rows, q = dim of the block, p fixed effects.
# Usage: REV_ARM=after|before Rscript dev/review-impsamp-perf.R G n q p N
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
a <- as.integer(commandArgs(TRUE))
G <- a[1]; n <- a[2]; q <- a[3]; p <- a[4]; N <- a[5]
ns <- asNamespace("frmtmb")
set.seed(1)
g <- factor(rep(seq_len(G), each = n))
X <- matrix(rnorm(G * n * p), G * n, p)
colnames(X) <- paste0("x", seq_len(p))
d <- data.frame(g = g, X)
eta <- 0.3 + 0.2 * X[, 1] + rnorm(G, 0, 0.5)[g]
if (q >= 2) eta <- eta + rnorm(G, 0, 0.3)[g] * X[, 1]
if (q >= 3) eta <- eta + rnorm(G, 0, 0.3)[g] * X[, 2]
d$y <- rpois(G * n, exp(eta))
re <- switch(q, "(1 | g)", "(1 + x1 | g)", "(1 + x1 + x2 | g)")
fo <- as.formula(paste("y ~", paste(colnames(X), collapse = " + "), "+", re))
fr <- frm(bf(fo), family = poisson(), data = d, dry_run = "frame")
nll <- ns$build_objective(fr)
tpl <- ns$make_start(fr, NULL, NULL)
lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
t_lapopt <- system.time(opt <- suppressWarnings(nlminb(lap$par, lap$fn, lap$gr,
  control = list(iter.max = 200))))[["elapsed"]]
at <- setNames(opt$par, names(lap$par))
lay <- ns$imp_layout(fr[["re_blocks"]])
gmap <- ns$imp_group_map(fr, lay)
otpl <- ns$imp_template(tpl, "b")
ctrl <- function() {
  t0 <- proc.time()[["elapsed"]]
  s <- 0
  for (i in 1:4e7) s <- s + i * 1e-9
  proc.time()[["elapsed"]] - t0
}
min_block <- function(f, rounds = 3L) {
  best <- Inf
  for (r in seq_len(rounds)) {
    k <- 0L
    t0 <- proc.time()[["elapsed"]]
    repeat {
      f()
      k <- k + 1L
      el <- proc.time()[["elapsed"]] - t0
      if (el > 1.2) break
    }
    best <- min(best, el / k)
  }
  best
}
c0 <- ctrl()
gc(reset = TRUE)
t_plan <- system.time(plan <- ns$imp_plan(lap, fr, lay, at, N, 1L))[["elapsed"]]
t_io <- system.time(io <- ns$build_importance_objective(fr, lay, gmap, plan))[["elapsed"]]
t_tape <- system.time(obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE))[["elapsed"]]
peak <- sum(gc()[, 6])
t_gr <- min_block(function() obj$gr(at))
t_fn <- min_block(function() obj$fn(at))
ad <- tryCatch(RTMB::MakeTape(function(x) 0, 1), error = function(e) NULL)
c1 <- ctrl()
cat(sprintf("RESULT\t%s\tG=%d n=%d q=%d p=%d N=%d\tplan=%.3f io=%.3f tape=%.3f gr=%.4f fn=%.4f peakMb=%.0f lapopt=%.2f control=%.2f,%.2f\n",
            arm, G, n, q, p, N, t_plan, t_io, t_tape, t_gr, t_fn, peak, t_lapopt, c0, c1))
if (length(a) >= 6 && a[6] == 1L) {
  w <- character(0)
  t_fit <- system.time(fit <- withCallingHandlers(
    tryCatch(frm(bf(fo), family = poisson(), data = d, importance = N),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }))[["elapsed"]]
  rr <- if (inherits(fit, "error")) "stopped" else fit$importance$rounds
  cat(sprintf("WHOLE\t%s\tG=%d n=%d q=%d p=%d N=%d\tfit=%.2f rounds=%s\n", arm, G,
              n, q, p, N, t_fit, rr))
}
