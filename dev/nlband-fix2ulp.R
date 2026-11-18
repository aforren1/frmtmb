# Lane nlband, fix round 2: the round's new assertions of
# test-fit-restart.R over 21 roundings (k = 0 to 20, a factor 1 + 2e-15 k):
#   SEP    separated bernoulli band data, seeds 2 to 4, x perturbed: no
#          hole message, a warning, the refusal message
#   STUDG  student() on the gaussian band data, seed 2, y perturbed: the
#          default fit has the hole message and a boundary message exactly
#          when the fit from (1, -2) has one
#   LOOSE  the short check's sweep case (n 60, y = 0.2 / (x - 2) + N(0, 2),
#          seed 13), y perturbed: the hole message, b in the hole, logLik
#          above the fit with nl_starts = FALSE
# Rscript dev/nlband-fix2ulp.R  (any R; the lane's library)
source("dev/nlband-common.R")
cat("R.home:", R.home(), "\n")
has <- function(r, cl) any(c(attr(r, "mc"), attr(r, "wc")) == cl) ||
  (cl == "frmtmb_boundary_fit" && any(attr(r, "mc") == "frmtmb_limit_fit"))
for (s in 2:4) for (k in 0:20) {
  d <- band_data(s)
  d$y <- as.integer(d$x > 2)
  d$x <- d$x * (1 + 2e-15 * k)
  r <- tt(frm(band_form, data = d, family = bernoulli()))
  ok <- !has(r, "frmtmb_nl_pole_start") && length(attr(r, "wc")) > 0 &&
    has(r, "frmtmb_nl_pole_refused")
  cat(sprintf("SEP seed %d k %d | %s | warn %s msg %s\n", s, k,
              if (ok) "PASS" else "FAIL", wsum(r), msum(r)))
}
d0 <- band_data(2)
for (k in 0:20) {
  d <- d0
  d$y <- d$y * (1 + 2e-15 * k)
  ref <- tt(frm(band_form, data = d, family = student(),
                start = list(beta = c(1, -2)),
                control = frmtmb_control(nl_starts = FALSE)))
  r <- tt(frm(band_form, data = d, family = student()))
  ok <- has(r, "frmtmb_nl_pole_start") &&
    has(r, "frmtmb_boundary_fit") == has(ref, "frmtmb_boundary_fit")
  cat(sprintf("STUDG k %d | %s | ref msg %s | msg %s\n", k,
              if (ok) "PASS" else "FAIL", msum(ref), msum(r)))
}
for (h in seq(0.2, 1.6, by = 0.02)) {
  set.seed(13)
  x <- runif(2000, 0, 4)
  x <- x[abs(x - 2) > h][seq_len(60)]
  hs <- ns$nl_pole_holes(sort(unique(-x)))
  if (!is.null(hs) && any(hs[, "lo"] < -2 & hs[, "hi"] > -2)) break
}
set.seed(513)
y0 <- 0.2 / (x - 2) + rnorm(60, 0, 2)
for (k in 0:20) {
  d <- data.frame(x = x, y = y0 * (1 + 2e-15 * k))
  off <- tt(frm(band_form, data = d,
                control = frmtmb_control(nl_starts = FALSE)))
  r <- tt(frm(band_form, data = d))
  b <- if (iserr(r)) NA else r$estimates$beta[[2L]]
  ok <- has(r, "frmtmb_nl_pole_start") &&
    isTRUE(b > -min(x[x > 2]) && b < -max(x[x < 2])) &&
    isTRUE(llv(r) > llv(off))
  cat(sprintf("LOOSE k %d | %s | ll %.6f off %.6f b %.4f msg %s\n", k,
              if (ok) "PASS" else "FAIL", llv(r), llv(off), b, msum(r)))
}
cat("DONE\n")
