# Load-independent instruments for test-perf.R's scaling test: tape
# nodes and bytes allocated while the objective is taped and evaluated
# once, at n and 100 n. Linear cost a + b n with a >= 0 bounds each
# ratio by 100. Also the instrument on a deliberately quadratic
# objective (an elementwise sub-assignment loop, dev/rtmb-pitfalls.md
# item 14), to see that it would catch one.
# Usage: Rscript dev/ciharden-perfcount.R <lib or base>
a <- commandArgs(TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
make <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
  d
}
alloc_bytes <- function(expr) {
  tf <- tempfile()
  Rprofmem(tf, threshold = 0)
  on.exit(Rprofmem(NULL))
  force(expr)
  Rprofmem(NULL)
  l <- readLines(tf)
  b <- suppressWarnings(as.numeric(sub(" *:.*", "", l)))
  sum(b, na.rm = TRUE)
}
census <- function(d) {
  t0 <- proc.time()[["elapsed"]]
  by <- alloc_bytes({
    fr <- frm(bf(y ~ x) + poisson(), data = d, dry_run = "frame")
    nll <- frmtmb:::build_objective(fr)
    obj <- RTMB::MakeADFun(nll, fr$par_template, silent = TRUE)
    obj$fn(obj$par)
    obj$gr(obj$par)
  })
  tp <- RTMB::GetTape(obj)
  nodes <- nrow(tp$data.frame())
  c(n = nrow(d), bytes = by, nodes = nodes,
    secs = proc.time()[["elapsed"]] - t0)
}
invisible(census(make(1000L, 70)))  # first call loads and caches
s <- census(make(1000L, 71))
l <- census(make(100000L, 72))
print(rbind(small = s, large = l))
cat(sprintf("ratio bytes %.2f nodes %.2f (bound 100)\n",
            l[["bytes"]] / s[["bytes"]], l[["nodes"]] / s[["nodes"]]))
# the same instrument on an objective with item 14's loop
quad <- function(n) {
  d <- make(n, 73)
  by <- alloc_bytes({
    f <- function(p) {
      "[<-" <- RTMB::ADoverload("[<-")
      eta <- p[1] + p[2] * d$x
      mu <- eta * 0
      for (i in seq_len(n)) mu[i] <- exp(eta[i])
      -sum(d$y * log(mu) - mu)
    }
    obj <- RTMB::MakeADFun(f, c(0, 0), silent = TRUE)
    obj$fn(obj$par)
  })
  by
}
q1 <- quad(1000L)
q2 <- quad(10000L)
cat(sprintf("loop objective: bytes at 1e3 %.4g, at 1e4 %.4g, ratio %.1f",
            q1, q2, q2 / q1), "for a tenfold n (linear bound 10)\n")
# the whole fit, as the wall-clock test timed it, with its evaluation
# count: allocation is setup + evaluations x per-evaluation cost
whole <- function(d) {
  by <- alloc_bytes(f <- suppressWarnings(frm(bf(y ~ x) + poisson(),
                                              data = d)))
  c(bytes = by, evals = f$opt$evaluations[["function"]],
    grads = f$opt$evaluations[["gradient"]])
}
invisible(whole(make(1000L, 70)))
ws <- whole(make(1000L, 71))
wl <- whole(make(100000L, 72))
print(rbind(small = ws, large = wl))
cat(sprintf("whole fit: ratio bytes %.2f\n", wl[["bytes"]] / ws[["bytes"]]))
