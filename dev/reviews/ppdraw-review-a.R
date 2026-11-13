# Review A: a zero inflation at its boundary with a predictor on zi
# (zi ~ x), and hurdle hu ~ x. Is the slope held? Does predict() agree with
# fitted()?
# Usage: Rscript ppdraw-review-a.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-ppdraw/dev/reviews/ppdraw-review-common.R")

mk <- function(zi, seed = 7, n = 200) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
  if (zi > 0) d$y[stats::runif(n) < zi] <- 0L
  d
}
show <- function(label, fit) {
  cat("\n==", label, "\n")
  fs <- flat_sides(fit)
  print(signif(fs, 4))
  r <- pz(fit)
  cat(sprintf("predict mean %.4f fitted mean %.4f z %.2f NA %.4f\n",
              r$pmean, r$fmean, r$z, r$na))
  cat("held message:", if (length(r$held)) r$held else "(none)", "\n")
  if (length(r$warns)) cat("warnings:", r$warns, sep = "\n  ")
  sdr <- frmtmb:::sdr_of(fit)
  cat("se_lost:", if (length(sdr$se_lost)) {
    paste(names(sdr$se_lost), sdr$se_lost, collapse = "; ")
  } else "(none)", "\n")
}

# A1: 3269's data with zi ~ x
d <- mk(0)
f1 <- q(frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = d))
show("A1 zib, 3269 data, zi ~ x", f1)

# A2: 3269's data, zi ~ x with x uncentered (x + 3)
d$x3 <- d$x + 3
f2 <- q(frm(bf(y | trials(N) ~ x, zi ~ x3), family = zero_inflated_binomial(),
            data = d))
show("A2 zib, 3269 data, zi ~ x + 3", f2)

# A3: zip on poisson data, zi ~ x
set.seed(11)
dp <- data.frame(x = stats::rnorm(300))
dp$y <- stats::rpois(300, exp(0.8 + 0.3 * dp$x))
f3 <- q(frm(bf(y ~ x, zi ~ x), family = zero_inflated_poisson(), data = dp))
show("A3 zip, poisson data, zi ~ x", f3)

# A4: hurdle_poisson with hu ~ x where no zeros at all
set.seed(12)
dh <- data.frame(x = stats::rnorm(300))
dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
cat("\nzeros in A4 data:", sum(dh$y == 0), "\n")
f4 <- q(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh))
show("A4 hurdle_poisson, no zeros, hu ~ x", f4)

# A5: zi ~ g, a factor on zi (two levels), 3269 data
d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
f5 <- q(frm(bf(y | trials(N) ~ x, zi ~ g), family = zero_inflated_binomial(),
            data = d))
show("A5 zib, 3269 data, zi ~ g", f5)
