# Lane nlband, fix round 1: the new assertions of test-fit-restart.R over
# 21 roundings (k = 0 to 20, a factor 1 + 2e-15 k): separated bernoulli
# band data, x perturbed (seeds 2 to 4; no hole message, a warning); and
# student band data, y perturbed (seeds 1 to 3; the hole message, no
# boundary message, logLik within grad_tol of the fit from (1, -2) on the
# same draw).
# Rscript dev/nlband-fix1ulp.R  (any R; the lane's library)
source("dev/nlband-common.R")
cat("R.home:", R.home(), "\n")
for (s in 2:4) for (k in 0:20) {
  d <- band_data(s)
  d$y <- as.integer(d$x > 2)
  d$x <- d$x * (1 + 2e-15 * k)
  r <- tt(frm(band_form, data = d, family = bernoulli()))
  ok <- !any(attr(r, "mc") == "frmtmb_nl_pole_start") &&
    length(attr(r, "wc")) > 0
  cat(sprintf("SEP seed %d k %d | %s | ll %.6g code %s warn %s msg %s\n", s, k,
              if (ok) "PASS" else "FAIL", llv(r),
              if (iserr(r)) "ERR" else r$opt$convergence, wsum(r), msum(r)))
}
for (s in 1:3) for (k in 0:20) {
  set.seed(s)
  x <- band_data(s)$x
  d <- data.frame(x = x, y = 1 / (x - 2) + 0.2 * rt(length(x), 3))
  d$y <- d$y * (1 + 2e-15 * k)
  ref <- tt(frm(band_form, data = d, family = student(),
                start = list(beta = c(1, -2)),
                control = frmtmb_control(nl_starts = FALSE)))
  r <- tt(frm(band_form, data = d, family = student()))
  cls <- c(attr(r, "mc"), attr(r, "wc"))
  ok <- any(cls == "frmtmb_nl_pole_start") &&
    !any(cls %in% c("frmtmb_limit_fit", "frmtmb_boundary_fit")) &&
    isTRUE(abs(llv(r) - llv(ref)) < 1e-3)
  cat(sprintf("STUD seed %d k %d | %s | ll %.8f ref %.8f warn %s msg %s\n", s,
              k, if (ok) "PASS" else "FAIL", llv(r), llv(ref), wsum(r),
              msum(r)))
}
cat("DONE\n")
