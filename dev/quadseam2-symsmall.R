# Fix round: SYMSMALL seed 1 of dev/reviews/quadseam2-review-esc2.R
# under one BLAS, with each escape_stationary() call logged (where,
# starts, gain, the alpha it ends at), to see where the fit ends.
#   Rscript dev/quadseam2-symsmall.R <libpath|base> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "BLAS:",
    extSoftVersion()[["BLAS"]], Sys.getenv("OPENBLAS_CORETYPE"), "\n")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1
ns <- asNamespace("frmtmb")
orig <- get("escape_stationary", ns)
unlockBinding("escape_stationary", ns)
assign("escape_stationary", function(obj, opt, frame, control, bounds,
                                     par_units = NULL, verbose = 0L,
                                     tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  st <- ns$stationary_escapes(obj, opt, frame)
  o <- orig(obj, opt, frame, control, bounds, par_units, verbose, tally)
  cat(sprintf("  ESC %s starts %d in %.6f out %.6f alpha_in %s alpha_out %s\n",
              where, length(st), opt$objective, o$objective,
              format(opt$par[names(opt$par) == "betad"][2], digits = 4),
              format(o$par[names(o$par) == "betad"][2], digits = 4)))
  o
}, envir = ns)
for (s in seeds) {
  set.seed(1000 + s)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:15, each = 4)))
  d$y <- 0.5 * d$x + rnorm(15, 0, 0.8)[d$g] + rnorm(60)
  w <- character()
  f <- withCallingHandlers(
    frm(bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ 1), family = skew_normal(),
        data = d, quadrature = TRUE, start = list(betad = c(0, 0))),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    })
  cat(sprintf("SEED %d ll %.10f alpha %s settle %s warn %s\n", s,
              as.numeric(logLik(f)), format(f$estimates$betad[2], digits = 6),
              paste(format(f$opt$quad_settle, digits = 3), collapse = "/"),
              paste(w, collapse = ",")))
}
