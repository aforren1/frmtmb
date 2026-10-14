# Lane illse: a covariate far from zero (dev/illse-designs.R des_offset,
# x on [1e5, 1e5 + 1], seed 2, no group term), on the base build and the
# worktree: code, Newton steps, log-likelihood against lm(), warnings.
#   Rscript dev/illse-offdbg.R base|dev [offset] [seed]
args <- commandArgs(TRUE)
arm <- args[1]
off <- if (length(args) > 1) as.numeric(args[2]) else 1e5
seed <- if (length(args) > 2) as.integer(args[3]) else 2L
.libPaths(c("C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-illse", quiet = TRUE))
} else suppressMessages(library(frmtmb))
set.seed(seed)
d <- data.frame(x = off + stats::runif(200))
d$y <- 1 + 2 * (d$x - off) + stats::rnorm(200, 0, 0.3)
w <- character()
f <- withCallingHandlers(frm(y ~ x, data = d, control = frmtmb_control(verbose = FALSE)),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
l <- lm(y ~ x, data = d)
cat("code", f$opt$convergence, f$opt$message, "| polish", f$opt$polish, "\n")
cat("logLik frm", as.numeric(logLik(f)), " lm", as.numeric(logLik(l)), "\n")
print(rbind(frm = fixef(f)[, 1], lm = coef(l)))
print(f$par_units)
cat(substr(w, 1, 300), sep = "\n")
