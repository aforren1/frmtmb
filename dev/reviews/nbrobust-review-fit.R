# Review of lane nbrobust, check 3: item 4903's data (seed 2, n = 200)
# fitted as negbinomial(), at a few roundings of x, from the default start
# and from log shape 30, under check_se "warning" and "ignore".
#   Rscript dev/reviews/nbrobust-review-fit.R <base|lib>
# Set OPENBLAS_CORETYPE / use the OpenBLAS R copy outside.
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
cat("La:", La_library(), " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
catch <- function(expr) {
  w <- character(0)
  m <- character(0)
  mc <- character(0)
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x))
    mc <<- c(mc, class(x)[1])
    invokeRestart("muffleMessage")
  })
  list(v = v, w = w, m = m, mc = mc)
}
set.seed(2)
d0 <- data.frame(x = stats::rnorm(200))
d0$y <- stats::rpois(200, exp(1 + 0.4 * d0$x))
for (k in c(0, 3, 7, 11, 15, 19)) {
  d <- d0
  d$x <- d$x * (1 + 2e-15 * k)
  for (st in c("default", "lsh30")) for (cse in c("warning", "ignore")) {
    ctl <- frmtmb_control(check_se = cse)
    r <- if (st == "default") {
      catch(frm(y ~ x, family = negbinomial(), data = d, control = ctl))
    } else {
      catch(frm(y ~ x, family = negbinomial(), data = d, control = ctl,
                start = list(betad = 30)))
    }
    f <- r$v
    gr <- tryCatch(max(abs(f$obj$gr(f$opt$par))), error = function(e) NA)
    cat(sprintf(paste("k %2d start %-7s cse %-7s code %d lsh %8.4f",
                      "logLik %.12f maxgr %9.3g\n"),
                k, st, cse, f$opt$convergence, f$opt$par[["betad"]],
                as.numeric(logLik(f)), gr))
    if (length(r$mc)) cat("   messages:", paste(unique(r$mc), collapse = ", "), "\n")
    for (w in unique(r$w)) cat("   WARNING:", substr(gsub("\n", " ", w), 1, 150), "\n")
  }
}
