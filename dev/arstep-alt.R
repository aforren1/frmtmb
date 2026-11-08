# Lane arstep: the second-start strategies on the step rows the seed scan
# found short (dev/arstep-seedscan.R), on the lane build. Strategy by the
# temporary option frmtmb.arstep_alt: none, local (the working marginal's
# next local maximum), scan (the objective scan without the working start),
# localscan (scan only where the marginal has a second maximum), both.
# Data as dev/review-ncgp-r3-ng.R (seed 74000 + 31 seed + n). Args: arm
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs(args[1])
f <- function(x) 2 * (x > 5) - 1
rows <- list(c("gar", 160, 1), c("gar", 160, 2), c("gar", 80, 22), c("gar", 80, 23), c("gar", 80, 4),
             c("student", 80, 1), c("student", 160, 1), c("student", 160, 22), c("student", 40, 5), c("student", 80, 6),
             c("student", 160, 4), c("gar", 160, 3), c("student", 80, 2))
for (rw in rows) {
  cs <- rw[1]; n <- as.integer(rw[2]); seed <- as.integer(rw[3])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- f(x)
  m <- switch(cs,
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
  for (how in c("none", "local", "scan", "localscan", "both")) {
    options(frmtmb.arstep_alt = how)
    t0 <- proc.time()[[3]]
    r <- rev_fit(m$fml, m$d, family = m$fam)
    el <- proc.time()[[3]] - t0
    if (inherits(r$fit, "error")) { cat(sprintf("ALT %s %d %d %s ERROR %s\n", cs, n, seed, how, conditionMessage(r$fit))); next }
    gs <- r$fit$opt$gp_second
    cat(sprintf("ALT %s %d %d %-9s ll %.4f conv %d w %d wall %.2f runs %s gain %s th %s\n", cs, n, seed, how, as.numeric(logLik(r$fit)),
                r$fit$opt$convergence, length(r$warn), el, if (is.null(gs)) "-" else gs[["runs"]],
                if (is.null(gs)) "-" else sprintf("%.4f", gs[["gain"]]), paste(signif(r$fit$estimates$theta, 4), collapse = ",")))
    flush.console()
  }
}
cat("DONE\n")
