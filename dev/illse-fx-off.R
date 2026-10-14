# Lane illse fix round 1: the offset fit of test-se-reparam.R
# (y ~ x, x on [1e5, 1e5 + 1]) under the Newton step changes.
#   Rscript dev/illse-fx-off.R base|lane|dev seed...
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
  d <- data.frame(x = 1e5 + stats::runif(200))
  d$y <- 1 + 2 * (d$x - 1e5) + stats::rnorm(200, 0, 0.3)
  w <- character()
  f <- withCallingHandlers(frm(y ~ x, data = d), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  ll <- as.numeric(logLik(lm(y ~ x, data = d)))
  cat(sprintf("seed %d | code %d | lm - fit %.3g | polish %s\n", seed,
              f$opt$convergence, ll - as.numeric(logLik(f)),
              paste(signif(f$opt$polish %||% NA, 4), collapse = ",")))
  cat(paste(" W:", substr(w, 1, 200)), sep = "\n")
}
