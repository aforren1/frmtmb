# Item 8: wall clock alone, with repeats. Args: arm(base|lane|warm) n rounds [family]
# Same design as dev/ncgp-perf.R: sin(x) + 0.3 e, seed 20261103.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; n <- as.integer(args[2]); rounds <- as.integer(args[3])
rev_libs(if (arm == "base") "base" else "lane")
source(file.path(REV_WT, "dev/review-ncgp-margin.R"))
control_block <- function() { m <- matrix(rnorm(1e6), 1000); t0 <- proc.time()[[3]]; for (i in 1:5) crossprod(m); (proc.time()[[3]] - t0) / 5 }
set.seed(20261103)
npos <- min(n, 500L)
pos <- sort(runif(npos, 0, 10))
x <- rep(pos, length.out = n)
d <- data.frame(x = x, y = sin(x) + 0.3 * rnorm(n))
for (i in seq_len(rounds)) {
  cb <- control_block()
  gc()
  st <- NULL; tm <- 0
  if (arm == "warm") {
    fr <- frm(y ~ gp(x), d, control = frmtmb_control(optCtrl = list(iter.max = 1L, eval.max = 2L), check_se = "ignore"))$frame
    t0 <- proc.time()[[3]]
    mo <- marg_opt(fr)
    tm <- proc.time()[[3]] - t0
    st <- mo$start
  }
  t0 <- proc.time()[[3]]
  fit <- suppressWarnings(if (is.null(st)) frm(y ~ gp(x), d) else frm(y ~ gp(x), d, start = st))
  el <- proc.time()[[3]] - t0
  t1 <- proc.time()[[3]]; s <- summary(fit); t2 <- proc.time()[[3]]
  cat(sprintf("RCO %s n %d pos %d round %d marginal %.2f frm %.2f total %.2f summary %.2f ll %.4f iter %s evals %s conv %d control %.3f\n",
              arm, n, npos, i, tm, el, tm + el, t2 - t1, logLik(fit), fit$opt$iterations,
              paste(fit$opt$evaluations, collapse = "/"), fit$opt$convergence, cb))
  flush.console()
}
cat("DONE\n")
