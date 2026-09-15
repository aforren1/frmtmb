# Lane `learnhier` (item 2.2): is the block the model fits actually
# non-Gaussian, and is the departure where the bias is?
#
#   Rscript dev/learnhier-gauss.R [n_seeds]
#
# NO FITS. Every number here comes from the drawn truths and the
# realized floors, which is the same construction the plan already
# relies on for the induced correlation. The data are re-simulated from
# the seed rather than read back, because the record stores the block's
# SUMMARY statistics and this needs the 100 deviations themselves.
#
# WHAT IS BEING TESTED, AND WHAT WOULD REFUTE IT. The rlddm arm's two
# failures are `sd(ndt)` biased low and `cor(bs, ndt)` driven toward -1,
# while `cor(alpha, ndt)` and `cor(drift, ndt)` are fine. The proposed
# explanation is that under `ndt_group(id)` the fitted deviation is
# `qlogis(ndt_i / floor_i)`, a logit of a ratio involving an OBSERVED
# MINIMUM, so the model is fitting the best Gaussian block to a
# non-Gaussian object. If the induced deviations come back near-Gaussian
# the explanation is wrong and the bias is something else.
#
# AND NON-GAUSSIANITY IS NECESSARY, NOT SUFFICIENT. If the deviations
# were OBSERVED, the Gaussian maximum likelihood estimate of the block
# would BE the sample moments, whatever their shape, so no departure
# could produce this bias. They are latent and integrated out, so two
# causes are in play at once: the shape of the block, and the ordinary
# downward bias of a maximum likelihood variance component estimated
# from limited per-group information. This script settles the first and
# cannot settle the second; the findings say so rather than claiming the
# whole mechanism.
#
# THE CONTROL IS THE POINT. Three of the four columns, `alpha`, `drift`
# and `bs`, are the drawn deviations untouched and are therefore
# Gaussian BY CONSTRUCTION. If the tests below flag them at more than
# their nominal rate, the instrument is broken and nothing it says about
# the fourth column can be believed.
source("dev/learnhier-env.R")
source("dev/learnhier-sim.R")
suppressPackageStartupMessages(library(frmtmb.learn))

a <- commandArgs(trailingOnly = TRUE)
NMAX <- if (length(a)) as.integer(a[[1L]]) else 1000L

fs <- list.files("dev/learnhier-rec", pattern = "^rlddm-[0-9]+[.]rds$",
                 full.names = TRUE)
seeds <- vapply(fs, function(f) readRDS(f)$seed, numeric(1))
seeds <- sort(unname(seeds))[seq_len(min(NMAX, length(fs)))]
cat("seeds ", length(seeds), ": ", min(seeds), " to ", max(seeds),
    "\n", sep = "")

skewness <- function(x) {
  z <- x - mean(x)
  mean(z^3) / mean(z^2)^1.5
}
ex_kurtosis <- function(x) {
  z <- x - mean(x)
  mean(z^4) / mean(z^2)^2 - 3
}

# Nonlinearity of the conditional mean of `y` given `x`: the p-value of
# the quadratic term. A bivariate Gaussian has a LINEAR conditional
# mean, so a quadratic that will not go away is a departure from
# bivariate normality that a covariance matrix cannot represent.
nonlin_p <- function(x, y) {
  m <- stats::lm(y ~ x + I(x^2))
  s <- summary(m)$coefficients
  if (nrow(s) < 3L) return(NA_real_)
  s[3L, 4L]
}

cols <- c("alpha", "drift", "bs", "ndt")
marg <- list()
pairs <- list()
for (sd_i in seeds) {
  d <- lh_rlddm_data(as.integer(sd_i))
  D <- attr(d, "dev_fitted")
  for (k in cols) {
    x <- D[, k]
    marg[[length(marg) + 1L]] <- data.frame(
      seed = sd_i, column = k, gaussian_by_construction = k != "ndt",
      skew = skewness(x), ex_kurt = ex_kurtosis(x),
      shapiro_p = stats::shapiro.test(x)$p.value,
      stringsAsFactors = FALSE)
  }
  pr <- lh_pairs(cols)
  for (j in seq_along(pr$lab)) {
    xa <- D[, pr$i[[j]]]
    xb <- D[, pr$j[[j]]]
    pairs[[length(pairs) + 1L]] <- data.frame(
      seed = sd_i, pair = pr$lab[[j]],
      involves_ndt = grepl("ndt", pr$lab[[j]], fixed = TRUE),
      pearson = stats::cor(xa, xb),
      spearman = stats::cor(xa, xb, method = "spearman"),
      nonlin_p = nonlin_p(xa, xb),
      stringsAsFactors = FALSE)
  }
  cat(".")
  utils::flush.console()
}
cat("\n")
marg <- do.call(rbind, marg)
pairs <- do.call(rbind, pairs)

cat("\n== marginal shape of the fitted-scale deviations, per column ==\n")
cat("alpha, drift and bs are the drawn deviations untouched and are the",
    " CONTROL.\n", sep = "")
mt <- do.call(rbind, lapply(cols, function(k) {
  s <- marg[marg$column == k, , drop = FALSE]
  data.frame(column = k, control = k != "ndt", n = nrow(s),
             mean_skew = mean(s$skew), mean_ex_kurt = mean(s$ex_kurt),
             shapiro_rejects_05 = sum(s$shapiro_p < 0.05),
             median_shapiro_p = stats::median(s$shapiro_p),
             stringsAsFactors = FALSE)
}))
mt[-(1:3)] <- lapply(mt[-(1:3)], function(x) signif(x, 4))
print(mt, row.names = FALSE)

cat("\n== pairwise: Pearson against Spearman, and a nonlinear",
    " conditional mean ==\n", sep = "")
pt <- do.call(rbind, lapply(unique(pairs$pair), function(p) {
  s <- pairs[pairs$pair == p, , drop = FALSE]
  data.frame(pair = p, involves_ndt = s$involves_ndt[[1L]], n = nrow(s),
             mean_pearson = mean(s$pearson),
             mean_spearman = mean(s$spearman),
             mean_gap = mean(abs(s$spearman) - abs(s$pearson)),
             nonlin_rejects_05 = sum(s$nonlin_p < 0.05, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
pt[-(1:3)] <- lapply(pt[-(1:3)], function(x) signif(x, 4))
print(pt, row.names = FALSE)

cat("\n== the mechanism check ==\n")
n <- length(seeds)
cat("Shapiro rejects at 0.05, control columns: ",
    sum(mt$shapiro_rejects_05[mt$control]), " of ", 3 * n,
    " (nominal would be about ", round(0.05 * 3 * n, 1), ")\n", sep = "")
cat("Shapiro rejects at 0.05, ndt column:      ",
    mt$shapiro_rejects_05[!mt$control], " of ", n, "\n", sep = "")
cat("nonlinear conditional mean, pairs WITHOUT ndt: ",
    sum(pt$nonlin_rejects_05[!pt$involves_ndt]), " of ", 3 * n, "\n",
    sep = "")
cat("nonlinear conditional mean, pairs WITH ndt:    ",
    sum(pt$nonlin_rejects_05[pt$involves_ndt]), " of ", 3 * n, "\n",
    sep = "")
saveRDS(list(marg = marg, pairs = pairs),
        "dev/learnhier-gauss.rds")
cat("\nwrote dev/learnhier-gauss.rds\n")
