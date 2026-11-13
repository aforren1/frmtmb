# Do the assertions of test-predict-draw-edge.R hold on every rounding?
# The covariate is perturbed by a few ulps, x * (1 + 2e-15 * k), k = 0..20,
# and each fit is asked what the test asks.
# Usage: Rscript dev/ppdraw-ulp.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
held_of <- function(expr) {
  said <- character(0)
  v <- withCallingHandlers(expr, frmtmb_draw_held = function(m) {
    said <<- c(said, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(v = v, said = said)
}
zib <- function(zi, k) {
  set.seed(7)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
  if (zi > 0) d$y[stats::runif(n) < zi] <- 0L
  d$x <- d$x * (1 + 2e-15 * k)
  d
}
for (zi in c(0, 0.3, 0.02)) {
  for (k in 0:20) {
    fit <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
               data = zib(zi, k))
    set.seed(1)
    r <- held_of(predict(fit, ndraws = 400, summary = FALSE))
    m <- rowMeans(r$v)
    z <- (mean(m) - mean(fitted(fit)[, "Estimate"])) /
      (stats::sd(m) / sqrt(400))
    cat(sprintf("zib zi=%.2f k=%2d zi_int=%9.4f held=%d z=%6.2f\n", zi, k,
                fit$estimates$betad[1], length(r$said), z))
  }
}
for (k in 0:20) {
  set.seed(7)
  d <- data.frame(x = stats::rnorm(400))
  kk <- stats::rbinom(400, 1, 0.3)
  d$yo <- cut(d$x + ifelse(kk == 1, 2, 0) + stats::rlogis(400),
              c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
  d$x <- d$x * (1 + 2e-15 * k)
  fit <- suppressWarnings(frm(bf(yo ~ x),
                              family = mixture(cumulative(), cumulative()),
                              data = d))
  set.seed(13)
  r <- held_of(predict(fit, newdata = data.frame(x = c(-1, 0, 1)),
                       summary = FALSE, ndraws = 4000))
  cat(sprintf("ordmix k=%2d tau2=%9.4f held=%d said_gap=%s NA=%d bad=%d\n",
              k, fit$estimates$tau_raw2[2], length(r$said),
              any(grepl("mu2_Intercept[2] - mu2_Intercept[1]", r$said,
                        fixed = TRUE)),
              sum(is.na(r$v)), sum(!r$v %in% 1:4)))
}
