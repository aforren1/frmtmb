# Item 2: can a finite-but-absurd predict() summary be told apart from a
# genuinely huge one? Candidate statistics per row, measured on fits
# where the huge value is GENUINE (heavy-tailed families, well
# identified) and on fits where a parameter sits in a nearly flat
# direction (an all-zero poisson cell, near-collinear predictors).
#   PREDFIX_ARM=base Rscript dev/predfix-absurd.R > dev/predfix-log/absurd-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
ND <- 1000L
stats_row <- function(dr, se_eta) {
  dr <- dr[is.finite(dr)]
  md <- stats::median(dr)
  a <- abs(dr - md)
  top <- sort(a, decreasing = TRUE)[seq_len(ceiling(length(a) / 100))]
  c(mean = mean(dr), sd = stats::sd(dr), median = md,
    # share of the total absolute deviation carried by the top 1% draws
    top1 = sum(top) / max(sum(a), .Machine$double.xmin),
    # mean against the 97.5% quantile: > 1 means a few draws drag the
    # mean past almost every draw
    mean_over_q975 = abs(mean(dr)) / max(abs(stats::quantile(dr, 0.975)),
                                         .Machine$double.xmin),
    se_eta = se_eta)
}
case <- function(label, kind, fit, nd) {
  set.seed(99)
  dr <- quiet(predict(fit, newdata = nd, ndraws = ND, summary = FALSE))
  se <- quiet(frm_linpred(fit, newdata = nd, se.fit = TRUE))$se.fit
  out <- lapply(seq_len(ncol(dr)), function(i) stats_row(dr[, i], se[i]))
  data.frame(label = label, kind = kind, row = seq_len(ncol(dr)),
             do.call(rbind, out))
}
rows <- list()
for (seed in 1:20) {
  set.seed(seed)
  n <- 200
  x <- rnorm(n)
  # GENUINE: well identified, heavy-tailed predictive law
  d <- data.frame(x = x, y = exp(rnorm(n, 1 + 0.3 * x, 2.5)))
  f <- quiet(frm(bf(y ~ x) + lognormal(), data = d))
  rows[[length(rows) + 1L]] <- case("lognormal sigma 2.5", "genuine", f,
                                    data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = 1 + 0.3 * x + rt(n, 1.3))
  f <- quiet(frm(bf(y ~ x) + student(), data = d))
  rows[[length(rows) + 1L]] <- case("student nu 1.3", "genuine", f,
                                    data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = rnbinom(n, mu = exp(2 + 0.3 * x), size = 0.15))
  f <- quiet(frm(bf(y ~ x) + negbinomial(), data = d))
  rows[[length(rows) + 1L]] <- case("negbinomial shape 0.15", "genuine", f,
                                    data.frame(x = c(0, 2)))
  # GENUINE extrapolation: the model's own answer far outside the data
  d <- data.frame(x = x, y = rpois(n, exp(0.3 + 0.4 * x)))
  f <- quiet(frm(bf(y ~ x) + poisson(), data = d))
  rows[[length(rows) + 1L]] <- case("poisson, x at 10 x the range",
                                    "extrapolation", f,
                                    data.frame(x = c(0, 30)))
  # FLAT: an all-zero cell of a poisson factor, rows the fit saw
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  yy <- rpois(n, exp(0.5 + 0.3 * x))
  yy[g == "c"] <- 0L
  d <- data.frame(x = x, g = g, y = yy)
  f <- quiet(frm(bf(y ~ x + g) + poisson(), data = d))
  rows[[length(rows) + 1L]] <- case("poisson, all-zero cell", "flat", f,
                                    data.frame(x = 0, g = c("a", "c")))
  # FLAT: near-collinear predictors, a row off their common direction
  x2 <- x + rnorm(n, 0, 1e-3)
  d <- data.frame(x = x, x2 = x2, y = rpois(n, exp(0.5 + 0.3 * x)))
  f <- quiet(frm(bf(y ~ x + x2) + poisson(), data = d))
  rows[[length(rows) + 1L]] <- case("poisson, x2 = x + N(0, 1e-3)", "flat",
                                    f, data.frame(x = c(0, 1), x2 = c(0, -1)))
}
res <- do.call(rbind, rows)
saveRDS(res, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/absurd.rds")
options(width = 200)
agg <- function(v, f) tapply(v, list(paste(res$label, "row", res$row)), f)
tab <- data.frame(
  kind = agg(res$kind, function(z) z[1]),
  n = agg(res$top1, length),
  est_max = agg(res$mean, max),
  esterr_max = agg(res$sd, max),
  top1_min = agg(res$top1, min), top1_max = agg(res$top1, max),
  mq_min = agg(res$mean_over_q975, min), mq_max = agg(res$mean_over_q975, max),
  se_eta_min = agg(res$se_eta, min), se_eta_max = agg(res$se_eta, max))
print(tab[order(tab$kind), ], digits = 4)
