# Reviewer of lane optrap, item 1: the unit is the family's start of
# sigma, which for gaussian() is sd(y), not the residual sd. When a
# covariate explains most of y the unit is far above the residual sd.
# y = 1 + k x + b_g + e (b sd 0.5, e sd 1), 20 groups of 5, seeds 1 to
# 40, k in 10, 100, 1e4; logLik against lme4 per build.
#   Rscript dev/optrap-rev-r2.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
X <- list()
for (k in c(10, 100, 1e4)) for (s in 1:40) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + k * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(100)
  w <- character()
  f <- withCallingHandlers(suppressMessages(frm(y ~ x + (1 | g), data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE,
    control = lmerControl(check.conv.singular = "ignore")))
  X[[length(X) + 1L]] <- data.frame(
    build = lab, k = k, seed = s, code = f$opt$convergence,
    dll = as.numeric(logLik(f)) - as.numeric(logLik(l4)),
    nwarn = length(w),
    unit = if (is.null(f$par_units)) 1 else f$par_units[1])
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-rev-log/r2-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
for (k in unique(X$k)) {
  x <- X[X$k == k, ]
  sh <- x$dll < -1e-4
  cat(sprintf(paste0("%s k %-6g unit %s | short of lme4 by > 1e-4: %d ",
                     "(min %.3g) | warned among short %d | code != 0 %d\n"),
              lab, k, paste(unique(x$unit), collapse = ","), sum(sh),
              min(x$dll), sum(x$nwarn[sh] > 0), sum(x$code != 0)))
}
