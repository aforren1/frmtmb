# Reviewer: recompute the reach-weighted rate and both intervals from
# results.tsv alone, with an independent estimator written from the
# design (PPS with replacement within strata).
r <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
fr <- read.delim("dev/equiv-data/strata.tsv")
cp <- read.delim("dev/equiv-frame-out/compat.tsv")
cp <- cp[cp$in_core_only, ]
w <- read.delim("dev/equiv-data/weights.tsv")
cp$w <- w$weight[match(cp$feature_a, w$feature)] *
  w$weight[match(cp$feature_b, w$feature)]
W <- tapply(cp$w, cp$status, sum) / sum(cp$w)
print(W)
stopifnot(sum(r$draws) == 200, nrow(r) == 172)
est <- function(ok) {
  # expand to draws; each draw is a Bernoulli trial within its stratum
  p <- v <- 0
  for (s in names(W)) {
    k <- r$status == s
    y <- rep(as.numeric(ok[k]), r$draws[k])
    n <- length(y)
    ph <- mean(y)
    p <- p + W[[s]] * ph
    # unbiased variance of a mean of iid draws
    v <- v + W[[s]]^2 * var(y) / n
  }
  c(est = 100 * p, lo = 100 * (p - 1.96 * sqrt(v)),
    hi = 100 * (p + 1.96 * sqrt(v)))
}
cat("agree, n-1 variance:\n"); print(round(est(r$agree), 2))
# the lane's form (divisor n)
p_n <- function(ok) {
  p <- v <- 0
  for (s in names(W)) {
    k <- r$status == s
    n <- sum(r$draws[k]); ph <- sum(r$draws[k] * ok[k]) / n
    p <- p + W[[s]] * ph; v <- v + W[[s]]^2 * ph * (1 - ph) / n
  }
  round(100 * c(p, p - 1.96 * sqrt(v), p + 1.96 * sqrt(v)), 2)
}
cat("agree, lane form:", p_n(r$agree), "\n")
cat("nodefect, lane form:", p_n(r$nodefect), "\n")
# nonparametric stratified bootstrap over draws
set.seed(1)
bs <- replicate(20000, {
  sum(vapply(names(W), function(s) {
    k <- which(r$status == s)
    y <- rep(as.numeric(r$agree[k]), r$draws[k])
    W[[s]] * mean(sample(y, length(y), replace = TRUE))
  }, 0))
})
cat("bootstrap (nonparametric):", round(100 * quantile(bs, c(.025, .975)), 2),
    "\n")
cat("\nverdict x stratum (draws):\n")
print(xtabs(draws ~ verdict + status, r))
cat("\nagree x verdict (unique):\n")
print(table(r$verdict, r$agree, useNA = "ifany"))
cat("\nagree, by top level (unique / draws):\n")
a <- r[r$agree, ]
print(table(a$top_level)); print(tapply(a$draws, a$top_level, sum))
cat("\nlevel string of agree pairs:\n"); print(table(a$level, useNA = "ifany"))
saveRDS(r, "dev/review-equiv-scratch/results.rds")
