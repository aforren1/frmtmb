# Reviewer of lane setier, re-check: the ri20 fits that stop short of
# lme4's maximum once the response is rescaled (dev/setier-rev2-scale2.R),
# and what each build tells the user about them.
#   Rscript dev/setier-rev2-short.R <lib> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:40
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
said <- function(w, m) {
  k <- c(if (any(grepl("^Optimizer did not report|gradient", w))) "CONV",
         if (any(grepl("^Standard errors are not available", w))) "SEWARN",
         if (any(grepl("^Boundary", m))) "BOUNDARYMSG")
  if (length(k)) paste(k, collapse = "+") else "nothing"
}
X <- list()
for (sc in c(1e-3, 1e3)) for (s in seeds) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * sc
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  dll <- as.numeric(logLik(f)) - as.numeric(logLik(l4))
  if (dll >= -1e-4) next
  X[[length(X) + 1L]] <- data.frame(scale = sc, seed = s,
                                    code = f$opt$convergence, dll = dll,
                                    lme4_singular = isSingular(l4),
                                    told = said(w, m))
}
X <- do.call(rbind, X)
print(X, row.names = FALSE, digits = 3)
print(table(scale = X$scale, told = X$told))
