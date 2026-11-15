# Review of lane wcheck2: the pair rule's cost on an UNWEIGHTED model
# with an observation-level effect and sigma ~ fac (no refusal possible:
# every excess is 0), against the whole frm() call.
#   Rscript dev/reviews/wcheck2-review-time2.R <base|lane>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
ns <- asNamespace("frmtmb")
for (cfg in list(c(1000, 50), c(2000, 100), c(4000, 200))) {
  nn <- cfg[1]; L <- cfg[2]
  set.seed(3)
  dd <- data.frame(obs = factor(seq_len(nn)), x = rnorm(nn),
                   fac = factor(sample(sprintf("l%04d", seq_len(L)), nn,
                                       replace = TRUE)))
  dd$y <- rnorm(nn)
  fl <- bf(y ~ x + (1 | obs), sigma ~ 1 + fac)
  fr <- q(frm(fl, data = dd, dry_run = "frame"))
  t0 <- proc.time()[["elapsed"]]
  r <- ns$check_weights_unbounded(fr$spec, fr, NULL)
  tc <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  f <- q(frm(fl, data = dd, control = frmtmb_control(check_olre = "ignore")))
  tf <- proc.time()[["elapsed"]] - t0
  cat(sprintf("unweighted n %d, %d levels: check %.2fs, whole frm() %.2fs\n",
              nn, L, tc, tf))
}
cat("DONE\n")
