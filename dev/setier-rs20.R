# Lane setier: what frmtmb says on dev/setier-singular.R's rs20 design
# (y ~ x + (1 + x | g), slope sd 0), seed given.
#   Rscript dev/setier-rs20.R <lib or "base"> <seed>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages({library(frmtmb); library(lme4)})
seed <- as.integer(args[2])
set.seed(seed)
d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
w <- character()
fit <- withCallingHandlers(frm(y ~ x + (1 + x | g), data = d),
  warning = function(c) {w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")},
  message = function(c) {w <<- c(w, paste("M:", conditionMessage(c)))
    invokeRestart("muffleMessage")})
cat("code", fit$opt$convergence, fit$opt$message, "evals", fit$opt$evals,
    "\n")
cat(substr(w, 1, 250), sep = "\n")
cat("par:", format(fit$opt$par, digits = 6), "\n")
cat("logLik frmtmb", as.numeric(logLik(fit)), "\n")
m <- lmer(y ~ x + (1 + x | g), data = d, REML = FALSE)
cat("logLik lme4", as.numeric(logLik(m)), " singular", isSingular(m), "\n")
print(VarCorr(m))
