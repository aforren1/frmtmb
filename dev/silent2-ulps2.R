# Lane silent2 fix round: the new assertions of test-ordinal-mixture.R
# ("a threshold with no standard error is said so" and "a degenerate fit
# with random effects names the starts to run") on 21 roundings of the
# data, x scaled by 1 + 2e-15 k, k = 0..20.
#   Rscript dev/silent2-ulps2.R
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
ws <- function(e) {
  w <- character()
  withCallingHandlers(e, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  w
}
cs1 <- function(k) {
  set.seed(1)
  n <- 300
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + rowSums(outer(lat, stats::quantile(lat, 1:3 / 4), ">"))
  d$x <- d$x * (1 + 2e-15 * k)
  d
}
flat <- function(seed, k) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d$x <- d$x * (1 + 2e-15 * k)
  d$g <- factor(rep_len(1:20, n))
  d
}
pa <- pb <- 0L
for (k in 0:20) {
  w <- ws(f <- frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()),
                   data = cs1(k)))
  fx <- suppressWarnings(fixef(f))
  lost <- any(is.nan(fx[grep("^mu2_Intercept", rownames(fx)), "Est.Error"]))
  dg <- grep("degenerate boundary", w, value = TRUE, fixed = TRUE)
  okA <- !(lost && any(grepl("its standard error is", dg, fixed = TRUE))) &&
    (!(lost && length(dg)) ||
       (any(grepl("has no standard error", dg, fixed = TRUE)) &&
          any(grepl("Standard errors are not available", w, fixed = TRUE))))
  w2 <- ws(frm(bf(y ~ x + z + (1 | g)), family = mixture(cumulative(),
                                                         sratio()),
               data = flat(9, k)))
  dg2 <- grep("degenerate boundary", w2, value = TRUE, fixed = TRUE)
  okB <- length(dg2) > 0L &&
    all(grepl("frmtmb_control(starts = 10)", dg2, fixed = TRUE))
  pa <- pa + okA
  pb <- pb + okB
  cat(sprintf("k %2d cs1 ll %.6f lost %s degen %d pass %s | re9 degen %d pass %s\n",
              k, as.numeric(logLik(f)), lost, length(dg), okA, length(dg2),
              okB))
}
cat(sprintf("SUMMARY cs1 pass %d of 21; re9 pass %d of 21\n", pa, pb))
