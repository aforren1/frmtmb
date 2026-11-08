# Cost of fitted() on new rows, per call, on one build (READERS_ARM). Each arm
# runs in its own process; the control is frm_linpred(), whose code this lane
# does not change, timed in the same process. Minimum of 5 rounds of blocks
# longer than 1.2 s. Data seeds 4101 and 9202 (readers-item1.R,
# readers-probe2.R).
source("dev/readers-setup.R")
blk <- function(f, reps) {
  t <- numeric(5)
  for (k in 1:5) t[k] <- system.time(for (i in seq_len(reps)) f())[["elapsed"]]
  min(t) / reps
}
set.seed(4101)
J <- 40; n <- 400
dd <- data.frame(x = runif(n, -1, 1), g = factor(sample(seq_len(J), n, TRUE)))
u <- rnorm(J)
dd$y <- rpois(n, exp(0.5 + 0.3 * dd$x + u[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd)
nd <- data.frame(x = c(-0.5, 0, 0.8), g = "new")
nd_seen <- transform(nd, g = dd$g[1])
cat(sprintf("control frm_linpred: %.2f ms\n",
            1000 * blk(function() frm_linpred(fit, newdata = nd_seen,
                                              type = "response"), 1500)))
cat(sprintf("fitted, seen level: %.1f ms\n",
            1000 * blk(function() fitted(fit, newdata = nd_seen), 40)))
cat(sprintf("fitted, new level: %.1f ms\n",
            1000 * blk(function() fitted(fit, newdata = nd,
                                         allow_new_levels = TRUE), 10)))
set.seed(9201)
d <- expand.grid(t = 1:25, s = factor(1:8))
d$x <- rnorm(nrow(d)); d$zo <- rnorm(nrow(d)); d$sz <- 0.5
e <- ave(rnorm(nrow(d)), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
zl <- d$zo + rnorm(nrow(d), 0, 0.5)
set.seed(9202)
lam <- exp(0.4 + 0.3 * d$x + 0.4 * zl + 0.4 * e)
cnt <- rpois(nrow(d), lam)
while (any(cnt == 0)) { i <- cnt == 0; cnt[i] <- rpois(sum(i), lam[i]) }
d$yh <- ifelse(runif(nrow(d)) < plogis(-0.5 + 0.8 * zl), 0L, cnt)
fh <- suppressWarnings(frm(bf(yh ~ x + me(zo, sz) + ar(time = t, gr = s,
                                                       cov = TRUE),
                              hu ~ me(zo, sz)), data = d,
                           family = hurdle_poisson()))
ndh <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1), zo = c(0, 1.5),
                  sz = 0.5)
cat(sprintf("control frm_linpred, hurdle: %.2f ms\n",
            1000 * blk(function() frm_linpred(fh, newdata = ndh,
                                              type = "response"), 1000)))
cat(sprintf("fitted, hurdle latent + me on two predictors: %.1f ms\n",
            1000 * blk(function() suppressWarnings(fitted(fh, newdata = ndh)),
                       3)))
# the reader marginaleffects' get_predict() calls, once per Jacobian column
ns <- asNamespace("frmtmb")
if (!is.null(ns$me_gh_moments)) {
  cat(sprintf("mean only, new level: %.1f ms\n", 1000 * blk(function() {
    ns$fitted_me_moments(fit, nd, NULL, NULL, NULL, TRUE, mean_only = TRUE)
  }, 10)))
  cat(sprintf("mean only, hurdle two directions: %.1f ms\n", 1000 * blk(
    function() ns$fitted_me_moments(fh, ndh, NULL, NULL, NULL, FALSE,
                                    mean_only = TRUE), 3)))
}
