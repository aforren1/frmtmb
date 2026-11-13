# Lane boundfam: dev/boundfam-cost.R's counters on fits where the
# screen asks (a negbinomial shape run off), with and without (1 | g).
#   Rscript dev/boundfam-cost2.R <lane lib>
a <- commandArgs(trailingOnly = TRUE)
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/boundfam-cost.R")
eval(parse(text = src[seq_len(grep("^cases <- list", src) - 1L)]))
cases <- list(
  nb_limit_fe = function() {
    set.seed(2)
    n <- 2e4
    d <- data.frame(x = rnorm(n))
    d$y <- rpois(n, exp(1 + 0.4 * d$x))
    frm(y ~ x, family = negbinomial(), data = d)
  },
  nb_limit_re = function() {
    set.seed(2)
    n <- 2e4
    d <- data.frame(x = rnorm(n), g = factor(sample(50, n, TRUE)))
    d$y <- rpois(n, exp(1 + 0.4 * d$x + rnorm(50, 0, 0.3)[d$g]))
    frm(y ~ x + (1 | g), family = negbinomial(), data = d)
  })
for (nm in names(cases)) {
  for (r in 1:3) {
    acc$t <- 0; acc$fn <- 0L; acc$gr <- 0L; acc$calls <- 0L
    t0 <- proc.time()[[3]]
    f <- suppressMessages(suppressWarnings(cases[[nm]]()))
    se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
    tt <- proc.time()[[3]] - t0
    cat(sprintf("%-12s rep %d: fit+vcov %.2f s, se_edge_limit %.2f s in %d calls, fn %d, gr %d, log shape %.3f, lost %s\n",
                nm, r, tt, acc$t, acc$calls, acc$fn, acc$gr,
                f$opt$par[["betad"]],
                paste(names(ns$sdr_of(f)$se_lost), collapse = ",")))
  }
}
