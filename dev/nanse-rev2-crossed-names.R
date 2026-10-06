# Reviewer, punch round 1: on the trial merge, how many coefficients
# does fixes' nl_flat_message() (with the projection rule) name on the
# crossed design of dev/nanse-rev-spread2.R, k = 20 and 60, and over
# 4 seeds? And is the SE warning then silent (one warning per fit)?
#   Rscript dev/nanse-rev2-crossed-names.R merge|release
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
for (k in c(20L, 60L)) {
  for (s in 3:6) {
    set.seed(s)
    d <- expand.grid(f = factor(seq_len(k)), g = factor(seq_len(k)))
    d <- d[sample(nrow(d), 6 * k), ]
    d$y <- rnorm(k)[d$f] + rnorm(k)[d$g] + rnorm(nrow(d), 0, 0.5)
    w <- character()
    fit <- withCallingHandlers(
      frm(bf(y ~ a + b, a ~ 0 + f, b ~ 0 + g, nl = TRUE), data = d),
      warning = function(x) {
        w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
      }, message = function(m) invokeRestart("muffleMessage"))
    fl <- grep("are not identified: at the optimum", w, value = TRUE)
    named <- if (length(fl)) {
      length(regmatches(fl, gregexpr("\\b(a_f|b_g)[0-9]+", fl))[[1]])
    } else 0L
    cat(sprintf(paste0("k = %d seed %d: warnings %d (flat %d, SE %d); ",
                       "coefficients named by the flat warning %d of %d; ",
                       "finite SEs %d\n"),
                k, s, length(w), length(fl),
                sum(grepl("Standard errors are not available", w)),
                named, 2L * k,
                sum(is.finite(suppressWarnings(fixef(fit)[, "Est.Error"])))))
  }
}
