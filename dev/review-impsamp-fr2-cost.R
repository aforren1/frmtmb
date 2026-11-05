# Fix round 2, item 7: what does the F2 search at the anchor cost? Time
# imp_prior_scales() against imp_plan() on large designs, by a block over 1.2 s.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
block <- function(f) {
  best <- Inf
  for (r in 1:3) {
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
for (cf in list(c(5000, 5, 2), c(5000, 5, 3), c(1000, 5, 1))) {
  G <- cf[1]; n <- cf[2]; q <- cf[3]
  set.seed(1)
  g <- factor(rep(seq_len(G), each = n))
  x1 <- rnorm(G * n); x2 <- rnorm(G * n)
  d <- data.frame(g, x1, x2, y = rpois(G * n, exp(0.3 + 0.2 * x1 + rnorm(G, 0, 0.5)[g])))
  re <- switch(q, "(1 | g)", "(1 + x1 | g)", "(1 + x1 + x2 | g)")
  fr <- frm(bf(as.formula(paste("y ~ x1 + x2 +", re))), family = poisson(), data = d, dry_run = "frame")
  nll <- ns$build_objective(fr)
  tpl <- ns$make_start(fr, NULL, NULL)
  lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
  opt <- suppressWarnings(nlminb(lap$par, lap$fn, lap$gr, control = list(iter.max = 200)))
  at <- setNames(opt$par, names(lap$par))
  lay <- ns$imp_layout(fr[["re_blocks"]])
  th <- unname(at[names(at) == "theta"])
  t_s <- block(function() ns$imp_prior_scales(lay, th))
  t_p <- system.time(plan <- ns$imp_plan(lap, fr, lay, at, 64L, 1L))[["elapsed"]]
  cat(sprintf("G %d q %d: imp_prior_scales %.3g s, imp_plan %.3f s, share %.3f percent\n", G, q, t_s, t_p, 100 * t_s / t_p))
}
