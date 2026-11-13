# Fix round, B1: what lane boundfam's build hands predict() for a "limit"
# parameter: its row of cov_fixed_prop (0, NaN?) and of the raw inverse.
# Usage: Rscript dev/ppdraw-fix-bf.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
show <- function(label, fit) {
  cat("\n==", label, "\n")
  sdr <- frmtmb:::sdr_of(fit)
  print(sdr$se_lost)
  cat("cov_fixed_prop:\n"); print(signif(sdr$cov_fixed_prop, 4))
  cat("cov_fixed_raw:\n"); print(signif(sdr$cov_fixed_raw, 4))
  cat("cov.fixed:\n"); print(signif(sdr$cov.fixed, 4))
  cat("root null:", is.null(sdr$cov_fixed_root), "\n")
}
set.seed(7)
d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
show("3269 zib", q(frm(bf(y | trials(N) ~ x),
                       family = zero_inflated_binomial(), data = d)))
show("3269 zib, zi ~ x", q(frm(bf(y | trials(N) ~ x, zi ~ x),
                               family = zero_inflated_binomial(), data = d)))
set.seed(7)
dg <- data.frame(x = stats::rnorm(200), s = stats::runif(200, 0.2, 1))
dg$y <- 1 + dg$x + stats::rnorm(200, 0, dg$s)
show("se(sigma = TRUE)", q(frm(bf(y | se(s, sigma = TRUE) ~ x),
                               family = gaussian(), data = dg)))
