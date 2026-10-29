## Lane rescorpp: the cost of the joint draw at 4000 draws x 500 rows x 3
## responses, against the per-response (independent) path.
##   Rscript dev/rescorpp-timing.R <lane|base>
## Arms run interleaved, 3 rounds, minimum reported; a fixed arithmetic
## control runs in each round to show the machine's load.
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-common.R")
set.seed(31); n <- 500
x <- rnorm(n)
R <- matrix(c(1, .8, .5, .8, 1, .3, .5, .3, 1), 3)
e <- MASS::mvrnorm(n, c(0, 0, 0), R)
d <- data.frame(x = x, y1 = 1 + x + e[, 1], y2 = -1 + 0.5 * x + e[, 2],
                y3 = 0.3 * x + e[, 3])
f <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE),
         family = gaussian(), data = d)
ds <- suppressWarnings(suppressMessages(frm_sample(
  f, chains = 2, iter = 4000, seed = 5, refresh = 0)))
cat("draws", ndraws(ds), " rows", n, "\n")
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
ctl_m <- matrix(rnorm(4e6), 2000)
arms <- list(
  control = function() invisible(crossprod(ctl_m)),
  independent = function() {
    for (r in c("y1", "y2", "y3")) posterior_predict(ds, resp = r)
  })
if (identical(arm, "lane")) {
  arms$joint <- function() posterior_predict(ds)
}
res <- matrix(NA_real_, 3, length(arms), dimnames = list(NULL, names(arms)))
for (k in 1:3) {
  for (a in names(arms)) {
    set.seed(k)
    res[k, a] <- tm(arms[[a]]())
  }
}
print(res)
cat("min seconds:\n"); print(apply(res, 2, min))
if (identical(arm, "lane")) {
  cat("joint / independent (min):",
      min(res[, "joint"]) / min(res[, "independent"]), "\n")
  ## where the joint path spends its time
  pf <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-timing.prof"
  Rprof(pf, interval = 0.005)
  invisible(posterior_predict(ds))
  Rprof(NULL)
  sp <- summaryRprof(pf)$by.total
  print(head(sp[, c("total.time", "total.pct")], 25))
  ## the draw alone, per call, on this size: one Cholesky per draw
  sh <- frmtmb.sample:::draws_fit_at(ds, 1L)
  dps <- frmtmb::eval_dpars(sh)
  ok <- rep(TRUE, n)
  t_draw <- tm(for (i in 1:4000) {
    frmtmb::predict_rescor_draw(sh, sh$spec$responses, dps, ok)
  })
  t_eval <- tm(for (i in 1:4000) frmtmb::eval_dpars(sh))
  # with the index the methods pass, as posterior_predict() calls it
  idx <- frmtmb.sample:::draws_index(ds)
  t_at <- tm(for (i in 1:4000) frmtmb.sample:::draws_fit_at(ds, i, idx))
  cat("per call ms: predict_rescor_draw", 1000 * t_draw / 4000,
      " eval_dpars", 1000 * t_eval / 4000,
      " draws_fit_at", 1000 * t_at / 4000, "\n")
}
