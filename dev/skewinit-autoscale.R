# The escape has to work with autoscale on, where the optimizer
# iterates in per-parameter units. A badly scaled covariate turns
# autoscale on.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
dd <- make_data(1, FALSE)
dd$xbig <- dd$xs * 1e4
m3 <- skew(dd$y)
bad <- list(betad = c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3))

one <- function(lab, fo, ...) {
  f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd, ...))
  e <- f$opt[["stationary_escape"]]
  cat(sprintf("%-34s ll %16.9f alpha %10.4f  units %s  esc %s\n", lab,
              as.numeric(logLik(f)),
              f$estimates$betad[["alpha_(Intercept)"]],
              if (is.null(f$par_units)) "none" else "on",
              if (is.null(e)) "-" else
                paste0(e[["starts"]], "/", format(e[["gain"]], digits = 5))))
}
fo_big <- bf(y ~ xbig, sigma ~ 1, alpha ~ 1)
one("xbig, default", fo_big)
one("xbig, raw-skew start", fo_big, start = bad)
one("xbig, autoscale TRUE", fo_big,
    control = frmtmb_control(autoscale = TRUE))
one("xbig, autoscale TRUE, raw start", fo_big, start = bad,
    control = frmtmb_control(autoscale = TRUE))
