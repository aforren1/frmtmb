# Short check, item 3, second part. Is theta's hold what makes predict() narrower?
# Theta free against theta held (mock of latres_held_pos), on fits whose theta is
# well determined. Coverage of a fresh count at each fitted row, given the true
# latent, and of a fresh series on new rows, by simulation from the truth.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
held_fn <- get("latres_held_pos", asNamespace("frmtmb"))
set_hold <- function(on) assignInNamespace("latres_held_pos", if (on) held_fn else function(fit, map) integer(0), ns = "frmtmb")
sim <- function(seed, s, ar = .4, G = 10, K = 6, b0 = .5) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- rnorm(n)
  d$ft <- factor(d$t)
  e <- as.vector(replicate(G, arima.sim(list(ar = ar), K, sd = s)))
  d$e <- e
  d$y <- rpois(n, exp(b0 + .3 * d$x + e))
  d
}
pred_both <- function(f, nd = NULL, anl = FALSE) {
  set_hold(TRUE); set.seed(1)
  a <- suppressMessages(predict(f, newdata = nd, ndraws = 800, allow_new_levels = anl))
  set_hold(FALSE); set.seed(1)
  b <- suppressMessages(predict(f, newdata = nd, ndraws = 800, allow_new_levels = anl))
  set_hold(TRUE)
  list(held = a, free = b)
}
cat("== Est.Error and width, held against free, ar1\n")
for (cc in list(c(11, .6), c(12, .6), c(21, 1), c(22, 1), c(7, .25), c(19, .25))) {
  d <- sim(cc[1], cc[2])
  f <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
  r <- frm(y ~ x + ar1(ft + 0 | g), family = poisson(), data = d)
  p <- pred_both(f)
  set.seed(1)
  pr <- suppressMessages(predict(r, ndraws = 800))
  vl <- sqrt(diag(as.matrix(vcov(f, full = TRUE))))
  cat(sprintf("seed %2d s %.2f se(log s) %6.2f: Est.Error held/ar1 %.3f free/ar1 %s held/free %s; max Estimate held %.3g free %.3g ar1 %.3g\n", cc[1], cc[2], vl["theta_1"],
              median(p$held[, 2] / pr[, 2]), format(median(p$free[, 2] / pr[, 2]), digits = 3), format(median(p$held[, 2] / p$free[, 2]), digits = 3),
              max(p$held[, 1]), max(p$free[, 1]), max(pr[, 1])))
}
cat("\n== coverage, 150 datasets per design: a fresh count at each fitted row given the true latent\n")
res <- list()
for (s in c(.25, .6, 1)) {
  hit <- matrix(NA, 150, 3)
  hitn <- matrix(NA, 150, 2)
  wid <- matrix(NA, 150, 3)
  for (i in 1:150) {
    d <- sim(6000 + i, s, G = 15, K = 8)
    f <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
    r <- tryCatch(frm(y ~ x + ar1(ft + 0 | g), family = poisson(), data = d), error = function(e) NULL)
    p <- pred_both(f)
    set.seed(6000 + i)
    y2 <- rpois(nrow(d), exp(.5 + .3 * d$x + d$e))
    hit[i, 1] <- mean(y2 >= p$held[, 3] & y2 <= p$held[, 4])
    hit[i, 2] <- mean(y2 >= p$free[, 3] & y2 <= p$free[, 4])
    wid[i, 1] <- median(p$held[, 4] - p$held[, 3])
    wid[i, 2] <- median(p$free[, 4] - p$free[, 3])
    if (!is.null(r)) {
      set.seed(1)
      pr <- suppressMessages(predict(r, ndraws = 800))
      hit[i, 3] <- mean(y2 >= pr[, 3] & y2 <= pr[, 4])
      wid[i, 3] <- median(pr[, 4] - pr[, 3])
    }
    nd <- data.frame(t = 1:8, g = factor(rep("new", 8)), x = 0, ft = factor(1:8, levels = levels(d$ft)))
    pn <- pred_both(f, nd, TRUE)
    yn <- rpois(8, exp(.5 + as.numeric(arima.sim(list(ar = .4), 8, sd = s))))
    hitn[i, 1] <- mean(yn >= pn$held[, 3] & yn <= pn$held[, 4])
    hitn[i, 2] <- mean(yn >= pn$free[, 3] & yn <= pn$free[, 4])
  }
  cat(sprintf("s %.2f fitted rows: coverage held %.3f free %.3f ar1 %.3f; median width held %.2f free %.2f ar1 %.2f | new series: held %.3f free %.3f\n",
              s, mean(hit[, 1]), mean(hit[, 2]), mean(hit[, 3], na.rm = TRUE), mean(wid[, 1]), mean(wid[, 2]), mean(wid[, 3], na.rm = TRUE), mean(hitn[, 1]), mean(hitn[, 2])))
}
