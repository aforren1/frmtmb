## Re-check, first priority: the shipped guard's false-alarm rate is
## now a number in the PLAN, "about 1 in 2,200", from a lognormal fit
## to 40 seeds with 0 of 40 at or below 1. That is an extrapolation
## into a tail no observation reaches. This attacks it three ways:
## does the shape hold, how far do OTHER shapes that fit the same body
## move the answer, and what does the data alone say without a shape.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ln <- readLines("dev/coh-paired.tsv", warn = FALSE)
ln <- ln[nzchar(trimws(ln))]
field <- function(k) {
  as.numeric(vapply(strsplit(ln, "\t", fixed = TRUE), function(p) {
    trimws(sub(paste0("^", k, "="), "",
               p[startsWith(p, paste0(k, "="))][1L]))
  }, character(1)))
}
q <- field("quotient")
n <- length(q)
lq <- log(q)
cat(sprintf("n = %d, quotient %.4f to %.4f, log mean %.4f sd %.4f\n",
            n, min(q), max(q), mean(lq), sd(lq)))

## 1. Does the shape hold? A normality test at n = 40 has little power,
##    so a pass licenses nothing; it is reported to show that it
##    licenses nothing, since the RAW scale passes too.
cat(sprintf("\nShapiro-Wilk on log(q): W = %.4f, p = %.3f\n",
            shapiro.test(lq)$statistic, shapiro.test(lq)$p.value))
cat(sprintf("Shapiro-Wilk on q:      W = %.4f, p = %.3f\n",
            shapiro.test(q)$statistic, shapiro.test(q)$p.value))
cat("Both fit the body, so the data do not choose between them.\n")

## 2. Shapes that fit the body, tails orders apart.
tail4 <- function(p) 1 - (1 - p)^4
p_lnorm <- pnorm(0, mean(lq), sd(lq))
p_norm <- pnorm(1, mean(q), sd(q))
sh <- mean(q)^2 / var(q)
p_gamma <- pgamma(1, shape = sh, rate = sh / mean(q))
## A t(5) on the log scale, matched on centre and scale: the same body,
## a tail that does not thin as fast.
sc <- sd(lq) / sqrt(5 / 3)
p_t5 <- pt((0 - mean(lq)) / sc, df = 5)
cat("\nfour shapes, all fitted to the same 40 points:\n")
cat(sprintf("%-22s per seed %-10.3g  min of four %-10.3g  1 in %.0f\n",
            c("lognormal (the lane's)", "normal on the raw scale",
              "gamma on the raw scale", "t(5) on the log scale"),
            c(p_lnorm, p_norm, p_gamma, p_t5),
            tail4(c(p_lnorm, p_norm, p_gamma, p_t5)),
            1 / tail4(c(p_lnorm, p_norm, p_gamma, p_t5))))
cat(sprintf("spread across shapes: a factor of %.0f\n",
            max(c(p_lnorm, p_norm, p_gamma, p_t5)) /
              min(c(p_lnorm, p_norm, p_gamma, p_t5))))

## 3. Parameter uncertainty inside the lane's OWN shape. The tail sits
##    3.7 sd out, so it is governed by an sd estimated from 40 points,
##    whose own standard error is sd / sqrt(2 (n - 1)).
cat(sprintf("\nse of the fitted sd: %.4f on an estimate of %.4f\n",
            sd(lq) / sqrt(2 * (n - 1)), sd(lq)))
set.seed(20260910)
B <- 20000L
bs <- replicate(B, {
  x <- rnorm(n, mean(lq), sd(lq))
  tail4(pnorm(0, mean(x), sd(x)))
})
cat(sprintf("lognormal refitted on 40 fresh draws: median 1 in %.0f, 95%% range 1 in %.0f to 1 in %.0f\n",
            1 / median(bs), 1 / quantile(bs, 0.975),
            1 / quantile(bs, 0.025)))
bn <- replicate(B, {
  x <- sample(lq, n, replace = TRUE)
  tail4(pnorm(0, mean(x), sd(x)))
})
cat(sprintf("nonparametric bootstrap of the same fit:  median 1 in %.0f, 95%% range 1 in %.0f to 1 in %.0f\n",
            1 / median(bn), 1 / quantile(bn, 0.975),
            1 / quantile(bn, 0.025)))

## 4. What the data say with no shape at all.
cp <- 1 - 0.05^(1 / n)
cat(sprintf("\ndistribution free, 0 of %d: Clopper-Pearson 95%% upper on the per-seed rate %.4f\n",
            n, cp))
cat(sprintf("  so the min of four fires at most %.3f of the time, about 1 in %.1f\n",
            tail4(cp), 1 / tail4(cp)))
cat(sprintf("  rule of three, 3/n: %.4f per seed\n", 3 / n))
