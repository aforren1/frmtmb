# Reviewer: tape cost of the summed CDFs (beta_binomial, compois) as
# the trials or the counts grow: time to build the objective and one
# gradient, and the tape size.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
set.seed(2)
n <- 200
for (tr in c(50, 500, 5000)) {
  d <- data.frame(x = rnorm(n), tr = tr)
  d$y <- rbinom(n, tr, rbeta(n, 4, 6))
  d$cc <- as.numeric(d$y > quantile(d$y, 0.8))
  t0 <- proc.time()[["elapsed"]]
  ob <- frm(y | trials(tr) + cens(cc) ~ x, family = beta_binomial(), data = d,
            dry_run = "objective")$obj
  t1 <- proc.time()[["elapsed"]]
  g <- ob$gr(ob$par)
  t2 <- proc.time()[["elapsed"]]
  d0 <- d; d0$cc <- 0
  ob0 <- frm(y | trials(tr) + cens(cc) ~ x, family = beta_binomial(), data = d0,
             dry_run = "objective")$obj
  t3 <- proc.time()[["elapsed"]]
  cat(sprintf("beta_binomial trials %5d, %d right-censored rows: build %.1f s, gradient %.2f s | uncensored build %.1f s\n",
              tr, sum(d$cc), t1 - t0, t2 - t1, t3 - t2))
}
for (lam in c(5, 50, 500)) {
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, lam)
  d$cc <- as.numeric(d$y > quantile(d$y, 0.8))
  t0 <- proc.time()[["elapsed"]]
  ob <- frm(y | cens(cc) ~ x, family = brmsfamily("com_poisson"), data = d,
            dry_run = "objective")$obj
  t1 <- proc.time()[["elapsed"]]
  g <- ob$gr(ob$par)
  t2 <- proc.time()[["elapsed"]]
  cat(sprintf("compois mean %4d: build %.1f s, gradient %.2f s\n", lam, t1 - t0, t2 - t1))
}
