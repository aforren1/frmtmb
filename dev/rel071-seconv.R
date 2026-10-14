# The case of test-se-conv.R's "an sd held near zero short of its
# maximum is not at the maximum" (lane illse, fix round 2) on one build:
# y ~ x + (1 | g) started at log sd -15 with iter.max = 3. Where does the
# stop end, with which code and which warnings, against the fit from the
# default start? It fails 5 of its expectations on the merged build.
#   Rscript dev/rel071-seconv.R <lib> [label]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("==", if (length(a) > 1) a[2] else "", "frmtmb",
    format(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
run <- function(expr) {
  ws <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = ws)
}
for (s in 1:4) {
  set.seed(s)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.4)
  ref <- run(frm(y ~ x + (1 | g), data = d))$v
  r <- run(frm(y ~ x + (1 | g), data = d, start = list(theta = -15),
               control = frmtmb_control(optCtrl = list(iter.max = 3,
                                                       eval.max = 1000))))
  f <- r$v
  cat(sprintf("seed %d: code %d (%s), iterations %s, logLik - default %.4f, theta %.3f\n",
              s, f$opt$convergence, f$opt$message, f$opt$iterations,
              as.numeric(logLik(f)) - as.numeric(logLik(ref)),
              f$estimates$theta[1]))
  if (!is.null(f$opt$restart)) {
    cat("   restart: code", f$opt$restart$convergence, f$opt$restart$message, "\n")
  }
  cat("   warnings:", length(r$w), "\n")
  for (w in r$w) cat("    -", substr(gsub("\n", " ", w), 1, 200), "\n")
  se <- tryCatch(sqrt(diag(vcov(f))), error = function(e) NA, warning = function(w) NA)
  cat("   SE(Intercept)", format(se[1], digits = 4), "default fit's",
      format(sqrt(diag(vcov(ref)))[1], digits = 4), "\n")
}
