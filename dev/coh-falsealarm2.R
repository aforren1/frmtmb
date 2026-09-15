## Punch round 2: what the 40 paired quotients can and cannot say about
## how often the shipped assertion would fire on correct data.
##
## The round 1 record quoted one shape's answer to two significant
## figures. This computes the model-free bound, which is what was
## OBSERVED, beside two fitted tails, so that the fit is visibly doing
## the work.
##
## Run: Rscript dev/coh-falsealarm2.R
ln <- readLines("dev/coh-paired.tsv", warn = FALSE)
ln <- ln[nzchar(trimws(ln))]
q <- as.numeric(vapply(strsplit(ln, "\t", fixed = TRUE), function(p) {
  v <- p[startsWith(p, "quotient=")]
  trimws(sub("^quotient=", "", v[1L]))
}, character(1)))
n <- length(q)
k <- 4L  ## the test takes the minimum of four seeds

cat(sprintf("n = %d, %d at or below 1\n", n, sum(q <= 1)))

## Model free. Clopper-Pearson is exact at 0 successes and is the
## bound to lead with, because it assumes nothing about the shape.
p_up <- 1 - 0.05^(1 / n)
cat(sprintf(paste0("Clopper-Pearson 95%% upper, per seed: %.4f",
                   "  (rule of three %.4f)\n"),
            p_up, 3 / n))
cat(sprintf("  so the four-seed statistic fires at most %.4f, 1 in %.1f\n",
            1 - (1 - p_up)^k, 1 / (1 - (1 - p_up)^k)))

## Two shapes the 40 points cannot choose between.
lq <- log(q)
cat(sprintf("\nShapiro-Wilk: log W %.4f p %.3f;  raw W %.4f p %.3f\n",
            stats::shapiro.test(lq)$statistic,
            stats::shapiro.test(lq)$p.value,
            stats::shapiro.test(q)$statistic,
            stats::shapiro.test(q)$p.value))
tail_of <- function(p1) 1 - (1 - p1)^k
p_lnorm <- stats::pnorm(0, mean(lq), stats::sd(lq))
p_norm <- stats::pnorm(1, mean(q), stats::sd(q))
cat(sprintf("lognormal on log: per seed %.3g, four-seed 1 in %.0f\n",
            p_lnorm, 1 / tail_of(p_lnorm)))
cat(sprintf("normal on raw:    per seed %.3g, four-seed 1 in %.0f\n",
            p_norm, 1 / tail_of(p_norm)))

## And inside the lane's own shape, the sd that governs the answer is
## itself estimated from 40 points.
set.seed(20260910)
B <- 20000L
r <- replicate(B, {
  s <- lq[sample.int(n, n, replace = TRUE)]
  1 / tail_of(stats::pnorm(0, mean(s), stats::sd(s)))
})
cat(sprintf(paste0("nonparametric bootstrap of the lognormal answer:",
                   " median 1 in %.0f, 95%% 1 in %.0f to 1 in %.0f\n"),
            stats::median(r), stats::quantile(r, 0.975),
            stats::quantile(r, 0.025)))
