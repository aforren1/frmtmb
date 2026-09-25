# Item 2, part 2: one candidate detector, calibrated. "one1" is the
# share of the row's total absolute deviation from the median carried by
# its single largest draw. The claim to test: when one1 > 1/2 the
# Estimate describes one draw, and it happens on flat directions and
# not on genuinely heavy tails. Also recorded: whether 0.62.0's
# non-finite warning already fires, and how far the Estimate moves
# between two predict() seeds (the direct measure of "absurd").
#   PREDFIX_ARM=base Rscript dev/predfix-absurd2.R > dev/predfix-log/absurd2-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
quiet <- function(expr) suppressMessages(expr)
ND <- 1000L
one_share <- function(v) {
  v <- v[is.finite(v)]
  a <- abs(v - stats::median(v))
  max(a) / max(sum(a), .Machine$double.xmin)
}
case <- function(label, kind, fit, nd) {
  warned <- FALSE
  run <- function(seed) {
    set.seed(seed)
    withCallingHandlers(
      quiet(predict(fit, newdata = nd, ndraws = ND, summary = FALSE)),
      warning = function(w) {
        if (grepl("not finite", conditionMessage(w))) warned <<- TRUE
        invokeRestart("muffleWarning")
      })
  }
  d1 <- run(99)
  d2 <- run(100)
  data.frame(label = label, kind = kind, row = seq_len(ncol(d1)),
             one1 = apply(d1, 2, one_share),
             est_ratio = apply(d1, 2, function(v) mean(v, na.rm = TRUE)) /
               apply(d2, 2, function(v) mean(v, na.rm = TRUE)),
             warned = warned)
}
rows <- list()
qfit <- function(...) suppressWarnings(quiet(frm(...)))
for (seed in 1:20) {
  set.seed(seed)
  n <- 200
  x <- rnorm(n)
  for (s in c(2.5, 3.5, 4.5)) {
    d <- data.frame(x = x, y = exp(rnorm(n, 1 + 0.3 * x, s)))
    rows[[length(rows) + 1L]] <- case(
      paste("lognormal sigma", s), "genuine",
      qfit(bf(y ~ x) + lognormal(), data = d), data.frame(x = c(0, 2)))
  }
  d <- data.frame(x = x, y = 1 + 0.3 * x + rt(n, 1.05))
  rows[[length(rows) + 1L]] <- case(
    "student, t(1.05) errors", "genuine",
    qfit(bf(y ~ x) + student(), data = d), data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = rnbinom(n, mu = exp(2 + 0.3 * x), size = 0.05))
  rows[[length(rows) + 1L]] <- case(
    "negbinomial shape 0.05", "genuine",
    qfit(bf(y ~ x) + negbinomial(), data = d), data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = rgamma(n, 0.2, 0.2 / exp(1 + 0.3 * x)))
  rows[[length(rows) + 1L]] <- case(
    "Gamma shape 0.2", "genuine",
    qfit(bf(y ~ x) + Gamma(link = "log"), data = d), data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = rpois(n, exp(0.3 + 0.4 * x)))
  fp <- qfit(bf(y ~ x) + poisson(), data = d)
  rows[[length(rows) + 1L]] <- case(
    "poisson, x at 10 and 100 x the range", "extrapolation", fp,
    data.frame(x = c(30, 300)))
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  yy <- rpois(n, exp(0.5 + 0.3 * x))
  yy[g == "c"] <- 0L
  rows[[length(rows) + 1L]] <- case(
    "poisson, all-zero cell", "flat",
    qfit(bf(y ~ x + g) + poisson(), data = data.frame(x = x, g = g, y = yy)),
    data.frame(x = 0, g = c("a", "c")))
  x2 <- x + rnorm(n, 0, 1e-3)
  rows[[length(rows) + 1L]] <- case(
    "poisson, x2 = x + N(0, 1e-3)", "flat",
    qfit(bf(y ~ x + x2) + poisson(),
         data = data.frame(x = x, x2 = x2, y = rpois(n, exp(0.5 + 0.3 * x)))),
    data.frame(x = c(0, 1), x2 = c(0, -1)))
  x3 <- x + rnorm(n, 0, 1e-2)
  rows[[length(rows) + 1L]] <- case(
    "poisson, x2 = x + N(0, 1e-2)", "flat",
    qfit(bf(y ~ x + x2) + poisson(),
         data = data.frame(x = x, x2 = x3, y = rpois(n, exp(0.5 + 0.3 * x)))),
    data.frame(x = c(0, 1), x2 = c(0, -1)))
}
res <- do.call(rbind, rows)
saveRDS(res, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/absurd2.rds")
options(width = 200)
key <- paste(res$label, "row", res$row)
tab <- data.frame(
  kind = tapply(res$kind, key, `[`, 1),
  n = tapply(res$one1, key, length),
  one1_min = tapply(res$one1, key, min),
  one1_med = tapply(res$one1, key, stats::median),
  one1_max = tapply(res$one1, key, max),
  over_half = tapply(res$one1 > 0.5, key, sum),
  est_ratio_min = tapply(res$est_ratio, key, min),
  est_ratio_max = tapply(res$est_ratio, key, max),
  warned = tapply(res$warned, key, sum))
print(tab[order(tab$kind), ], digits = 4)
