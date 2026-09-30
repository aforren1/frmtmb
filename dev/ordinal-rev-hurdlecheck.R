# Reviewer: brms's hurdle_cumulative nthres and top-category branch, on the
# worker's data (seed 20260930) and the reviewer's (seed 20261001).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yh <- ifelse(runif(n) < plogis(-1 + 0.3 * d$x), 0L, d$y)
cat("worker yh table:", table(d$yh), "\n")
for (fam in list(brms::hurdle_cumulative(threshold = "sum_to_zero"),
                 brms::hurdle_cumulative("probit", threshold = "equidistant"),
                 brms::hurdle_cumulative())) {
  sd <- brms::standata(brms::bf(yh ~ x), data = d, family = fam)
  sc <- brms::stancode(brms::bf(yh ~ x), data = d, family = fam)
  L <- strsplit(sc, "\n")[[1]]
  cat(fam$link, fam$threshold, ": nthres", sd$nthres, "| top branch:",
      trimws(grep("else if [(]y == nthres", L, value = TRUE)[1]), "\n")
}
for (f in list(brms::bf(yh ~ x, hu ~ x), brms::bf(yh ~ x), brms::bf(yh ~ x, disc ~ 0 + z))) {
  sc <- brms::stancode(f, data = d, family = brms::hurdle_cumulative(threshold = "sum_to_zero"))
  L <- strsplit(sc, "\n")[[1]]
  cat(deparse(f$formula), ":", trimws(grep("target [+]=", L, value = TRUE)[1]), "\n")
  i <- grep("real hurdle_cumulative_logit_lpmf", L)
  cat(L[i:(i + 16)], sep = "\n")
}
