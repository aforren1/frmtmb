# Lane optrap, consolidation fixes: lane seflags' B5 design without the
# group term (degree-8 raw powers of x on [1, 2] plus z and w, seeds 4
# and 5), where the 0.70.0 merge lost every SE, z, wb and sigma among
# them, because the Newton steps ended where the convergence check
# skips the SE check. SEs of z, wb and log sigma against the same model
# with poly(x, 8), the fit's warnings, and its Newton steps.
#   Rscript dev/optrap-c3-b5.R <lib> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else c(4, 5)
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "\n")
for (s in seeds) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200),
                  w = factor(sample(c("a", "b"), 200, TRUE)))
  d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
    stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  w <- character()
  f <- withCallingHandlers(
    frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + w, data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  fr <- suppressWarnings(frm(y ~ poly(x, 8) + z + w, data = d))
  sa <- suppressWarnings(summary(f)$coefficients)
  sb <- summary(fr)$coefficients
  rat <- c(z = sa$mu["z", "Std. Error"] / sb$mu["z", "Std. Error"],
           wb = sa$mu["wb", "Std. Error"] / sb$mu["wb", "Std. Error"],
           sigma = sa$sigma[1, "Std. Error"] / sb$sigma[1, "Std. Error"])
  cat(sprintf("seed %d code %d | logLik - poly %.3g | steps %s | SE ratio z %.4f wb %.4f sigma %.4f | warnings %d: %s\n",
              s, f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(fr)),
              (f$opt$polish %||% c(steps = 0))[["steps"]], rat[1], rat[2],
              rat[3], length(w), substr(paste(w, collapse = " | "), 1, 160)))
}
