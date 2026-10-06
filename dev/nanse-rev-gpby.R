# Reviewer: the two fits of test-gp-by.R "check C" that the trial merge
# warns on. Are their SEs lost on the release build too?
#   Rscript dev/nanse-rev-gpby.R release|merge
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
src <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat"
for (h in list.files(src, "^helper.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = globalenv())
}
gpby_data <- function(n = 60, seed = 5) {
  set.seed(seed)
  d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                  f = factor(rep(c("a", "b", "c"), length.out = n)),
                  w = stats::runif(n, 0.5, 2))
  d$y <- 0.5 + ifelse(d$f == "a", sin(d$x),
                      ifelse(d$f == "b", cos(d$x), 0.2 * d$x)) +
    stats::rnorm(n, 0, 0.3)
  d
}
d <- gpby_data()
set.seed(31)
d$y2 <- 1 + cos(d$x) * (d$f == "a") + stats::rnorm(nrow(d), 0, 0.4)
d$ysig <- 0.5 + sin(d$x) +
  stats::rnorm(nrow(d), 0, exp(-1 + 0.3 * (d$f == "b") * sin(d$x)))
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(v = v, w = w)
}
fits <- list(
  fs = function() frm(bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)), data = d),
  fb = function() frm(mvbf(bf(y ~ gp(x, by = f, k = 8)),
                           bf(y2 ~ gp(x, by = f, k = 6))) +
                        set_rescor(FALSE), data = d, family = gaussian()))
for (nm in names(fits)) {
  r <- cap(fits[[nm]]())
  f <- r$v
  sdr <- RTMB::sdreport(f$obj)
  se0 <- sqrt(pmax(diag(sdr$cov.fixed), NA))
  nmz <- frmtmb:::outer_par_names(f)
  cat("\n==", nm, "code", f$opt$convergence, "\n")
  cat("  plain sdreport SE:", paste(sprintf("%s=%.3g", nmz, se0),
                                    collapse = " "), "\n")
  cat("  estimates:", paste(sprintf("%s=%.3g", nmz, f$opt$par),
                            collapse = " "), "\n")
  for (x in r$w) cat("  warn:", substr(x, 1, 200), "\n")
  H <- optimHess(f$opt$par, f$obj$fn, f$obj$gr)
  cat("  |H| row max:", paste(sprintf("%s=%.2g", nmz,
                                      apply(abs(H), 1, max)),
                              collapse = " "), "\n")
}
