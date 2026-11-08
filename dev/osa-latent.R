# Lane osa: defect 3, bernoulli residuals(type = "osa") near -3.7 on latent
# fits. The construction of dev/review-latres-osa2.R (seeds 1 to 3, 15
# groups of 10, latent AR(1) at 0.5 with sd 0.6), run on this lane's
# libraries.
#
#   Rscript dev/osa-latent.R before|after
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " frmtmb from", dirname(find.package("frmtmb")), "\n")
sumr <- function(r) {
  sprintf("mean %6.3f sd %5.3f min %6.2f max %6.2f nan %d ks.p %.3g",
          mean(r, na.rm = TRUE), sd(r, na.rm = TRUE), min(r, na.rm = TRUE),
          max(r, na.rm = TRUE), sum(!is.finite(r)),
          suppressWarnings(stats::ks.test(r[is.finite(r)], "pnorm")$p.value))
}
cat("BEGIN osa-latent", arm, "\n")
for (seed in 1:3) {
  set.seed(seed)
  d <- expand.grid(t = 1:10, g = factor(1:15))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(15, arima.sim(list(ar = .5), 10, sd = .6)))
  d$id <- factor(seq_len(n))
  d$ft <- factor(d$t)
  d$bin <- rbinom(n, 1, plogis(.2 + .4 * d$x + e))
  d$bi <- rbinom(n, 5, plogis(.2 + .4 * d$x + e))
  d$tr <- 5L
  d$cnt <- rpois(n, exp(.5 + .4 * d$x + e))
  d$nb <- rnbinom(n, size = 3, mu = exp(.5 + .4 * d$x + e))
  d$ord <- factor(findInterval(.3 * d$x + e + rlogis(n), c(-1, 0, 1)) + 1L,
                  ordered = TRUE)
  specs <- list(
    bern_latent = list(bin ~ x + ar(t, g, cov = TRUE), bernoulli()),
    bern_re = list(bin ~ x + (1 | g), bernoulli()),
    bern_ar1re = list(bin ~ x + ar1(ft + 0 | g), bernoulli()),
    bern_plain = list(bin ~ x, bernoulli()),
    binom_latent = list(bi | trials(tr) ~ x + ar(t, g, cov = TRUE),
                        binomial()),
    binom_re = list(bi | trials(tr) ~ x + (1 | g), binomial()),
    pois_latent = list(cnt ~ x + ar(t, g, cov = TRUE), poisson()),
    nb_latent = list(nb ~ x + ar(t, g, cov = TRUE), negbinomial()),
    ord_latent = list(ord ~ x + ar(t, g, cov = TRUE), cumulative()),
    ord_re = list(ord ~ x + (1 | g), cumulative()))
  for (nm in names(specs)) {
    f <- suppressWarnings(frm(specs[[nm]][[1]], family = specs[[nm]][[2]],
                              data = d))
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(suppressWarnings(
      as.matrix(residuals(f, type = "osa"))[, 1]), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("seed %d %-13s %s  %.2fs\n", seed, nm,
                if (inherits(r, "error")) {
                  paste("ERR", substr(conditionMessage(r), 1, 100))
                } else sumr(r), el))
  }
}
cat("END osa-latent\n")
