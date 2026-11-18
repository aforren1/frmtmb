# Lane ordstarts, build round: the new assertions of test-mixture-starts.R
# on every rounding, x scaled by 1 + 2e-15 k, k = 0..20 (lane-env-r084.md).
# 1. "the ranking keeps today's three proper fits proper": cum_sratio
#    n = 1000 seeds 13 and 34, n = 300 seed 28, starts = 10.
# 2. "a mixture of two families keeps a swapped fit that is better":
#    cum_sratio n = 1000 seed 13, the default.
#   Rscript dev/ordstarts-ulps.R <which: rank|swap> [ks]
a <- commandArgs(trailingOnly = TRUE)
which_one <- a[1]
ks <- if (length(a) > 1) eval(parse(text = a[2])) else 0:20
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
oms_data <- function(seed, n = 300, k = 0) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d$x <- d$x * (1 + 2e-15 * k)
  d
}
fit <- function(d, ctl) {
  w <- character()
  m <- character()
  f <- withCallingHandlers(
    frm(bf(y ~ x + z), family = mixture(cumulative(), sratio()), data = d,
        control = ctl),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, class(x)[1])
      invokeRestart("muffleMessage")
    })
  list(f = f, w = w, m = m)
}
ns <- asNamespace("frmtmb")
pass <- 0L
tot <- 0L
if (which_one == "rank") {
  for (cfg in list(c(1000, 13), c(1000, 34), c(300, 28))) {
    for (k in ks) {
      d <- oms_data(cfg[2], cfg[1], k)
      off <- fit(d, frmtmb_control(starts = 0))
      on <- fit(d, frmtmb_control(starts = 10))
      cl <- ns$mixture_ord_start_class(on$f, "y")
      ok <- !any(grepl("degenerate boundary", off$w)) &&
        !any(grepl("degenerate boundary", on$w)) &&
        as.numeric(logLik(on$f)) >= as.numeric(logLik(off$f)) - 1e-3 &&
        !length(cl$step) && !length(cl$boundary)
      tot <- tot + 1L
      pass <- pass + ok
      cat(sprintf("n %d seed %d k %2d off %.4f on %.4f ok %s\n", cfg[1],
                  cfg[2], k, as.numeric(logLik(off$f)),
                  as.numeric(logLik(on$f)), ok))
    }
  }
} else {
  for (k in ks) {
    d <- oms_data(13, 1000, k)
    off <- fit(d, frmtmb_control(starts = 0))
    on <- fit(d, frmtmb_control())
    ms <- on$f$opt$mixture_starts
    llo <- as.numeric(logLik(off$f))
    lln <- as.numeric(logLik(on$f))
    ok <- ms[["swaps"]] == 1 && ms[["starts"]] == 0 && lln >= llo - 1e-3
    if (ms[["swapped"]] == 1) {
      cl <- ns$mixture_ord_start_class(on$f, "y")
      ok <- ok && sum(on$m == "frmtmb_mixture_swap") == 1L &&
        lln > llo + 1e-3 && on$f$opt$convergence == 0L &&
        !length(cl$step) && !length(cl$boundary) &&
        !any(grepl("degenerate boundary", on$w))
    } else {
      ok <- ok && !any(on$m == "frmtmb_mixture_swap") &&
        identical(on$f$opt$par, off$f$opt$par)
    }
    tot <- tot + 1L
    pass <- pass + ok
    cat(sprintf("k %2d off %.4f on %.4f swapped %d ok %s\n", k, llo, lln,
                ms[["swapped"]], ok))
  }
}
cat("SUMMARY", which_one, "pass", pass, "of", tot, "\n")
