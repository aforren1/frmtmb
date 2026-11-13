# Review of lane ordstart: the three-component cases of the sweep where
# the lane replaced a step function by a "flat boundary" point. What is
# that point: converged? where does the optimizer go when continued from
# it? what does mixture_ord_degeneracy() measure there?
#   Rscript ordstart-review-cum3.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/reviews/ordstart-review-sweep.R")
eval(parse(text = src[grep("^gen <- ", src):(grep("^cases <- list", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
fam <- mixture(cumulative(), cumulative(), cumulative())
op <- options(width = 160)
for (s in c(1, 4, 6, 7)) {
  d <- gen(s, n = 300, K = 5)
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x + z), family = fam, data = d),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  f0 <- q(frm(bf(y ~ x + z), family = fam, data = d,
              control = frmtmb_control(starts = 0)))
  cat(sprintf(paste("\nSEED %d lane ll %.4f code %d (%s) iter %s max|par| %.4g",
                    "| starts=0 ll %.4f code %d max|par| %.4g\n"), s,
              as.numeric(logLik(f)), f$opt$convergence, f$opt$message,
              f$opt$iterations %||% NA, max(abs(f$opt$par)),
              as.numeric(logLik(f0)), f0$opt$convergence,
              max(abs(f0$opt$par))))
  cat("ms:", paste(names(f$opt$mixture_starts), f$opt$mixture_starts,
                   collapse = " "), "\n")
  print(frmtmb:::mixture_ord_degeneracy(f, "y"), digits = 4)
  cat("largest |par| of the lane's point:\n")
  p <- f$opt$par
  print(signif(p[order(-abs(p))][1:6], 4))
  # continue the optimization from the lane's point, no starts
  g <- q(frm(bf(y ~ x + z), family = fam, data = d, start = f$estimates,
             control = frmtmb_control(starts = 0)))
  dg <- frmtmb:::mixture_ord_degeneracy(g, "y")
  cat(sprintf("CONTINUED from the lane's point: ll %.4f code %d max|par| %.4g max sharpen %s\n",
              as.numeric(logLik(g)), g$opt$convergence, max(abs(g$opt$par)),
              paste(signif(dg$sharpen, 3), collapse = ",")))
}
cat("DONE\n")
