# Reviewer re-check: the cost of the laplace probe and watch on a large
# posterior_predict(). 20000 rows in 200 groups, 1000 draws (data seed 5,
# draws seed 1). The full-draws object runs neither the probe nor the
# watch, and the two objects hold the same outer values.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-rev-r2-01-watch.R")
eval(parse(text = src[seq_len(grep("set.seed(77)", src, fixed = TRUE)[1L] - 1L)]))
set.seed(5)
G <- 200; n <- 100
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(G, 0, 0.6)[dd$g] + rnorm(nrow(dd), 0, 0.7)
fit <- suppressMessages(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
p <- lap_pair(fit, n = 1000L)
tm <- function(expr) { gc(); t <- system.time(expr)[["elapsed"]]; t }
reps <- 3L
out <- matrix(NA_real_, reps, 4, dimnames = list(NULL,
  c("predict_full", "predict_lap", "epred_full", "epred_lap")))
for (r in seq_len(reps)) {
  set.seed(1); out[r, 1] <- tm(a <- posterior_predict(p$full, re_formula = NA))
  set.seed(1); out[r, 2] <- tm(b <- posterior_predict(p$lap, re_formula = NA))
  out[r, 3] <- tm(ea <- posterior_epred(p$full, re_formula = NA))
  out[r, 4] <- tm(eb <- posterior_epred(p$lap, re_formula = NA))
}
cat("identical predict:", identical(a, b), " identical epred:", identical(ea, eb), "\n")
print(out)
med <- apply(out, 2, median)
cat(sprintf("median seconds: predict full %.2f, laplace %.2f (ratio %.3f); epred full %.2f, laplace %.2f (ratio %.3f)\n",
            med[1], med[2], med[2] / med[1], med[3], med[4], med[4] / med[3]))
# the watch alone, on dpars the size posterior_predict() hands it
w <- frmtmb.sample:::draws_laplace_watch(p$lap, "posterior_predict()")
dp <- list(mu = rnorm(nrow(dd)), sigma = rep(0.7, nrow(dd)))
tw <- system.time(for (i in 1:1000) w(dp))[["elapsed"]]
cat(sprintf("watch alone, 1000 calls on 2 x %d cells: %.3f s (%.2f ms per draw)\n",
            nrow(dd), tw, 1000 * tw / 1000))
