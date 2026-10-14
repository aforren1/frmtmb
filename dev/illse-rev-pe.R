# Reviewer of lane illse, attack 4 (rtmb-pitfalls item 21): the new
# route's results across processes pinned to P and E cores. Fits the
# lane's own test designs and saves what a test reads.
#   Rscript dev/illse-rev-pe.R <out.rds>
Sys.setenv(REV_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
out <- commandArgs(TRUE)[1]
res <- list()
mkp <- function(seed, re) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, length.out = 200)))
  eta <- sin(3 * d$x)
  if (re) eta <- eta + stats::rnorm(20, 0, 0.5)[d$g]
  d$y <- eta + stats::rnorm(200, 0, 0.3)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  d
}
fo8 <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
nd <- data.frame(x = c(1, 1.5, 2), g = factor(1, levels = 1:20))
for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
for (s in 1:3) {
  d <- mkp(s, FALSE)
  f <- suppressWarnings(frm(fo8, data = d))
  res[[paste("deg8", s)]] <- list(par = f$opt$par, code = f$opt$convergence,
    polish = f$opt$polish, V = suppressWarnings(vcov(f, full = TRUE)),
    fit = fitted(f, newdata = nd, re_formula = NA)[, "Est.Error"])
  d <- mkp(s, TRUE)
  f <- suppressWarnings(frm(update(fo8, . ~ . + (1 | g)), data = d))
  res[[paste("deg8re", s)]] <- list(par = f$opt$par, code = f$opt$convergence,
    polish = f$opt$polish, V = suppressWarnings(vcov(f, full = TRUE)),
    fit = fitted(f, newdata = nd)[, "Est.Error"])
}
set.seed(2)
d <- data.frame(x = 1e5 + stats::runif(200))
d$y <- 1 + 2 * (d$x - 1e5) + stats::rnorm(200, 0, 0.3)
f <- suppressWarnings(frm(y ~ x, data = d))
res[["offset"]] <- list(par = f$opt$par, code = f$opt$convergence,
  polish = f$opt$polish, V = suppressWarnings(vcov(f, full = TRUE)))
saveRDS(res, out)
cat("done\n")
