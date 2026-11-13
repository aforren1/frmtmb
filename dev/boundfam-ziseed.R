# Lane boundfam fix round: seeds where a zero_inflated_poisson's zi
# intercept runs to its limit on plain poisson data (a quiet control).
.libPaths(c("C:/Users/adf44/source/r/wt-boundfam-lib", "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (s in 1:12) {
  set.seed(s)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(1 + 0.3 * d$x))
  w <- 0L; m <- character()
  f <- withCallingHandlers(frm(bf(y ~ x, zi ~ 1), family = zero_inflated_poisson(), data = d),
    warning = function(e) { w <<- w + 1L; invokeRestart("muffleWarning") },
    message = function(e) { m <<- c(m, class(e)[1]); invokeRestart("muffleMessage") })
  cat(sprintf("seed %d: zi %.3f lost %s warnings %d messages %s\n", s,
              f$opt$par[3], paste(ns$sdr_of(f)$se_lost, collapse = ","), w,
              paste(m, collapse = ",")))
}
