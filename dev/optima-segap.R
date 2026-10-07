# Lane optima: replacement fixtures for the SE-check tests that read
# the softmax plateau of mo(), which the lane removes. A monotonic
# predictor with an unobserved interior category: its two steps enter
# only as their sum, so the simplex has an exactly flat direction on
# every data set, not on a knife edge. What does the SE check report?
#   Rscript dev/optima-segap.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
se_mo_gap <- function(seed) {
  set.seed(seed)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
  data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
}
for (fo in list(ls ~ mo(inc), ls ~ mo(inc) + age, ls ~ mo(inc) * age)) {
  for (s in 1:3) {
    w <- character()
    f <- withCallingHandlers(frm(fo, data = se_mo_gap(s)),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
    sdr <- frmtmb:::sdr_of(f)
    cat(deparse(fo), "seed", s, ": code", f$opt$convergence, "; warnings",
        length(w), "; lost", paste(names(sdr$se_lost), sdr$se_lost,
                                   collapse = ", "), "\n")
    if (length(w)) cat("   ", substr(w, 1, 200), "\n")
  }
}
# D = 2 with the middle category unobserved: the one coordinate does
# not enter the likelihood at all
for (s in 1:5) {
  set.seed(s)
  d2 <- data.frame(inc = sample(c(0L, 2L), 100, TRUE), age = rnorm(100, 40, 10))
  d2$ls <- c(30, 50, 70)[d2$inc + 1L] + rnorm(100, sd = 7)
  for (fo in list(ls ~ mo(inc), ls ~ mo(inc) * age)) {
    w <- character()
    f <- withCallingHandlers(frm(fo, data = d2), warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
    H <- f$obj$he(f$opt$par)
    z <- grep("^zeta", names(f$opt$par))
    cat(deparse(fo), "seed", s, ": warnings", length(w), "; lost",
        paste(names(frmtmb:::sdr_of(f)$se_lost), collapse = ","),
        "; max |H row of zeta|", format(max(abs(H[z, ])), digits = 3),
        "; gradient", format(f$obj$gr(f$opt$par)[z], digits = 3), "\n")
  }
}
