# Lane silent2: the new ordinal-mixture assertions of test-ordinal-mixture.R
# on 21 roundings of the data (x scaled by 1 + 2e-15 k, k = 0..20), as
# lane-env asks. Per k: cum_sratio seed 16 must warn exactly when its
# largest threshold standard error (fixef()) is above 100; seed 2 under
# both families must not warn, and its SE must stay below 100.
#   Rscript dev/silent2-ulps.R base|lane
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
gen <- function(seed, k = 0, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d$x <- d$x * (1 + 2e-15 * k)
  d
}
one <- function(d, fam) {
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x + z), family = fam, data = d),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  fx <- suppressWarnings(fixef(f))
  se <- max(fx[grep("_Intercept[[]", rownames(fx)), "Est.Error"])
  list(degen = any(grepl("degenerate boundary", w, fixed = TRUE)), se = se,
       ll = as.numeric(logLik(f)))
}
ok16 <- ok2s <- ok2c <- 0L
for (k in 0:20) {
  a <- one(gen(16, k), mixture(cumulative(), sratio()))
  b <- one(gen(2, k), mixture(cumulative(), sratio()))
  c <- one(gen(2, k), mixture(cumulative(), cumulative()))
  p16 <- identical(a$degen, a$se > 100)
  ok16 <- ok16 + p16
  ok2s <- ok2s + (!b$degen && b$se < 100)
  ok2c <- ok2c + (!c$degen && c$se < 100)
  cat(sprintf("k %2d seed16 degen %s se %.4g ll %.6f pass %s | seed2 sr %s %.3g | cc %s %.3g\n",
              k, a$degen, a$se, a$ll, p16, b$degen, b$se, c$degen, c$se))
}
cat(sprintf("SUMMARY seed16 pass %d of 21; seed2 sratio pass %d of 21; seed2 cumulative pass %d of 21\n",
            ok16, ok2s, ok2c))
