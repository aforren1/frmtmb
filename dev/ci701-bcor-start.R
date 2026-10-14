# Lane ci701: the boundary-correlation fit of test-review-v29.R
# (set.seed(61)) from moved starts, to find natural fits on the side the
# macOS runner took (the optimizer reports convergence and the
# standard-error check runs). Start: theta = (0, 0, t0) with t0 on a
# grid; prints the code, whether a convergence warning came, the lost
# set, the boundary message and the sd rows of confint_varcorr().
#
#   Rscript dev/ci701-bcor-start.R <lib|base>
args <- commandArgs(TRUE)
lib <- args[1]
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
for (t0 in c(0.5, 1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100, 200, 500, 1000,
             3000)) {
  w <- character(0); bm <- 0L
  fit <- withCallingHandlers(
    frm(f, data = dd, start = list(theta = c(0, 0, t0))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    frmtmb_boundary_fit = function(m) {
      bm <<- bm + 1L; invokeRestart("muffleMessage")
    })
  lost <- frmtmb:::sdr_of(fit)$se_lost
  ci <- suppressWarnings(confint_varcorr(fit))
  s <- ci[ci$type == "sd", ]
  cat(sprintf(
    "RES t0=%g code=%d t=%.1f nwarn=%d bmsg=%d lost=%s sd_lwr=%s sd_upr=%s\n",
    t0, fit$opt$convergence, fit$estimates$theta[3], length(w), bm,
    if (length(lost)) paste(names(lost), lost, sep = ":", collapse = ",")
    else "none",
    paste(sprintf("%.4f", s$lwr), collapse = ","),
    paste(sprintf("%.4f", s$upr), collapse = ",")))
}
