# Lane illse fix round 2, R2-1: the reviewer's real case of
# dev/illse-rev2-flat2.R (y ~ x + (1 | g), true sd 0.5, start log sd
# -15, iter.max 3): the convergence warnings in full.
#   Rscript dev/illse-fx-r2flat.R base|lane|dev seed...
args <- commandArgs(TRUE)
arm <- args[1]
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-illse-lib"),
            file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
for (seed in as.integer(args[-1])) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.4)
  w <- character()
  f <- withCallingHandlers(
    frm(y ~ x + (1 | g), data = d, start = list(theta = -15),
        control = frmtmb_control(optCtrl = list(iter.max = 3,
                                                eval.max = 1000))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  ref <- suppressWarnings(frm(y ~ x + (1 | g), data = d))
  cat(sprintf("seed %d | code %d | short by %.4g\n", seed,
              f$opt$convergence,
              as.numeric(logLik(ref)) - as.numeric(logLik(f))))
  cat(paste(" W:", substr(w, 1, 230)), sep = "\n")
}
