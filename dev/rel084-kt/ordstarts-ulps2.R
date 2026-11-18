# rel084-kt copy of dev/ordstarts-ulps2.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Lane ordstarts, build fix round 1: the new assertions of
# test-mixture-starts.R on every rounding, x scaled by 1 + 2e-15 k,
# k = 0..20 (lane-env-r084.md), as dev/ordstarts-ulps.R does.
# weights: "constant weights do not change the maximum returned"
#   (cum_sratio n = 1000 seed 13, weights 1/1000, 1, 10: the same swap
#   verdict, and the logLik per unit weight within 1e-3 of the gain the
#   swap buys at weights 1 over the fit with starts = 0; no parameter
#   check, as in the test since the short check of 2026-11-18).
# seflat: "the starts named by the warning run on a default flat by its
#   SE" (cum_cum n = 300 seeds 3 and 16, starts = 10: the starts ran and
#   the fit is not below the default's).
#   Rscript dev/ordstarts-ulps2.R <weights|seflat> [ks]
a <- commandArgs(trailingOnly = TRUE)
which_one <- a[1]
ks <- if (length(a) > 1) eval(parse(text = a[2])) else 0:20
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "; BLAS",
    extSoftVersion()[["BLAS"]], "\n")
oms_data <- function(seed, n = 300, k = 0) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d$x <- d$x * (1 + 2e-15 * k)
  d
}
q <- function(e) suppressWarnings(suppressMessages(e))
pass <- 0L
tot <- 0L
if (which_one == "weights") {
  for (k in ks) {
    d <- oms_data(13, 1000, k)
    got <- lapply(c(1 / 1000, 1, 10), function(w) {
      d$w <- w
      f <- q(frm(bf(y | weights(w) ~ x + z),
                 family = mixture(cumulative(), sratio()), data = d))
      list(ll = as.numeric(logLik(f)) / w,
           swapped = f$opt$mixture_starts[["swapped"]])
    })
    d$w <- 1
    off <- q(frm(bf(y | weights(w) ~ x + z),
                 family = mixture(cumulative(), sratio()), data = d,
                 control = frmtmb_control(starts = 0)))
    gain <- got[[2L]]$ll - as.numeric(logLik(off))
    ok <- gain > 0 && all(vapply(got[-2L], function(g) {
      identical(g$swapped, got[[2L]]$swapped) &&
        abs(g$ll - got[[2L]]$ll) < 1e-3 * gain
    }, NA))
    tot <- tot + 1L
    pass <- pass + ok
    cat(sprintf("k %2d ll %.6f %.6f %.6f swapped %s ok %s\n", k,
                got[[1]]$ll, got[[2]]$ll, got[[3]]$ll,
                paste(vapply(got, `[[`, 0, "swapped"), collapse = ""), ok))
  }
} else {
  for (s in c(3, 16)) {
    for (k in ks) {
      d <- oms_data(s, 300, k)
      fam <- mixture(cumulative(), cumulative())
      off <- q(frm(bf(y ~ x + z), family = fam, data = d))
      on <- q(frm(bf(y ~ x + z), family = fam, data = d,
                  control = frmtmb_control(starts = 10)))
      ok <- identical(on$opt$mixture_starts[["starts"]], 10) &&
        as.numeric(logLik(on)) >= as.numeric(logLik(off)) - 1e-3
      tot <- tot + 1L
      pass <- pass + ok
      cat(sprintf("seed %d k %2d off %.4f on %.4f ok %s\n", s, k,
                  as.numeric(logLik(off)), as.numeric(logLik(on)), ok))
    }
  }
}
cat("SUMMARY", which_one, "pass", pass, "of", tot, "\n")
