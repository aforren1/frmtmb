# Lane silent2 fix round, 5214 (2): what building the covariance at the
# fit costs for an ordinal mixture with random effects whose
# standard-error check would wait. The reviewer's design
# (dev/ordstart-cost.R re1000): n = 1000, 50 groups, y ~ x + z + (1 | g),
# mixture(cumulative(), cumulative()), seeds 1 to 3. Per seed and round:
# the default fit's elapsed seconds (on the build before the forcing, the
# check defers), then se_fit_analysis(force = TRUE) on a copy of that
# fit, which is what the forcing adds. Minimum over 3 rounds of each.
#   Rscript dev/silent2-recost.R [seeds]
args <- commandArgs(trailingOnly = TRUE)
seeds <- if (length(args)) eval(parse(text = args[1])) else 1:3
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
mk <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
ns <- asNamespace("frmtmb")
for (s in seeds) {
  d <- mk(s, n = 1000)
  d$g <- factor(rep_len(1:50, 1000))
  fit_t <- an_t <- numeric(0)
  for (r in 1:3) {
    t0 <- proc.time()[["elapsed"]]
    f <- quiet(frm(bf(y ~ x + z + (1 | g)),
                   family = mixture(cumulative(), cumulative()), data = d))
    fit_t <- c(fit_t, proc.time()[["elapsed"]] - t0)
    deferred <- is.null(f$cache$se_fit_an) && is.null(f$cache$sdr)
    t0 <- proc.time()[["elapsed"]]
    an <- quiet(ns$se_fit_analysis(f, force = TRUE))
    an_t <- c(an_t, proc.time()[["elapsed"]] - t0)
  }
  cat(sprintf(paste("RECOST seed %d np %d evals %s deferred %s fit %.2fs",
                    "covariance %.2fs ratio (fit + cov) / fit %.3f\n"),
              s, length(f$opt$par), format(f$opt$evals %||% NA), deferred,
              min(fit_t), min(an_t), (min(fit_t) + min(an_t)) / min(fit_t)))
}
