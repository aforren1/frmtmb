# Lane nbrobust: a negbinomial fit on poisson data started at a large log
# shape, under check_se = "ignore" (no SE analysis, so no limit rule to
# read the gradient), over ulp perturbations of the data.
#   Rscript dev/nbrobust-start.R <base|lib> [ls0] [cse]
# Data: seed 2, n = 200 (backlog item 4903). x * (1 + 2e-15 * k),
# k = 0..20.
args <- commandArgs(TRUE)
arm <- args[1]
ls0 <- if (length(args) >= 2) as.numeric(args[2]) else 30
cse <- if (length(args) >= 3) args[3] else "ignore"
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), " R:", R.home(),
    " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
conds <- function(expr) {
  w <- list()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
set.seed(2)
d0 <- data.frame(x = stats::rnorm(200))
d0$y <- stats::rpois(200, exp(1 + 0.4 * d0$x))
for (k in 0:20) {
  d <- d0
  d$x <- d$x * (1 + 2e-15 * k)
  r <- conds(frm(y ~ x, family = negbinomial(), data = d,
                 start = list(betad = ls0),
                 control = frmtmb_control(check_se = cse)))
  fit <- r$value
  p <- fit$opt$par
  wm <- vapply(r$warnings, conditionMessage, "")
  cat(sprintf(paste("RESULT start k %2d code %d (%s) shape %.4f logLik",
                    "%.12f gradwarn %s nwarn %d nmsg %d\n"), k,
              fit$opt$convergence, fit$opt$message, p[length(p)],
              as.numeric(logLik(fit)),
              any(grepl("Large maximum absolute gradient", wm)),
              length(wm), length(r$messages)))
  for (w in wm) cat("   W", substr(w, 1, 200), "\n")
}
cat("END\n")
