# Reviewer of lane ordstarts, build round. Probe 5 (question 3): does a
# fit create .Random.seed where none existed? One fresh process per
# library; each arm starts with no .Random.seed.
#   Rscript ordstarts-build-review-p5.R <base|lane>
which_lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (which_lib == "lane") "C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
d <- mk(13, 1000)
ge <- globalenv()
arm <- function(label, expr) {
  if (exists(".Random.seed", envir = ge)) rm(".Random.seed", envir = ge)
  v <- suppressWarnings(suppressMessages(expr))
  ms <- v$opt$mixture_starts
  cat(sprintf("%-40s .Random.seed exists after: %s | record: %s\n", label,
              exists(".Random.seed", envir = ge),
              if (is.null(ms)) "NULL" else
                paste(names(ms), signif(ms, 3), sep = "=", collapse = " ")))
}
arm("gaussian lm-like", frm(bf(x ~ z), family = gaussian(), data = d))
arm("cumulative alone", frm(bf(y ~ x + z), family = cumulative(), data = d))
arm("cum_sratio starts = 0", frm(bf(y ~ x + z),
                                 family = mixture(cumulative(), sratio()),
                                 data = d,
                                 control = frmtmb_control(starts = 0)))
arm("cum_sratio default", frm(bf(y ~ x + z),
                              family = mixture(cumulative(), sratio()),
                              data = d))
arm("cum_sratio starts = 10", frm(bf(y ~ x + z),
                                  family = mixture(cumulative(), sratio()),
                                  data = d,
                                  control = frmtmb_control(starts = 10)))
cat("DONE p5\n")
