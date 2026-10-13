# Reviewer of lane optrap, re-check: lane moparity's own frame for
# mo(x) * f (its library, read only), with optrap's predictor_col_spread()
# (sourced from this tree) applied to its mu predictor.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
ex <- parse("R/autoscale.R")
env <- new.env(parent = asNamespace("frmtmb"))
for (e in ex) if (is.call(e) && identical(deparse(e[[2]]), "predictor_col_spread"))
  eval(e, env)
set.seed(1)
n <- 120
lev <- c("a", "b", "c", "d", "e")
d <- data.frame(x = factor(sample(lev, n, TRUE), levels = lev, ordered = TRUE),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$y <- as.integer(d$x) * (1 + (d$f == "b")) + rnorm(n)
for (fo in list(y ~ mo(x) * f, y ~ mo(x):f, y ~ mo(x) + (mo(x) | f))) {
  fr <- tryCatch(frm(fo, data = d, dry_run = "frame"), error = function(e) e)
  if (inherits(fr, "error")) { cat(deparse(fo), "ERROR", conditionMessage(fr), "\n"); next }
  lp <- Filter(function(l) identical(l$dpar, "mu"), fr$linpreds)[[1]]
  cat(deparse(fo), "| columns", colnames(lp$X), "\n")
  for (m in lp$mo) cat("   col", m$col, "codes", length(m$codes), "mult",
                       if (is.null(m$mult)) "NULL" else length(m$mult), "\n")
  cat("   spreads", env$predictor_col_spread(lp), "\n")
}
d$f[3] <- NA
fr <- tryCatch(frm(y ~ mo(x) * f, data = d, dry_run = "frame"), error = function(e) e)
if (inherits(fr, "error")) cat("NA in f: ERROR", conditionMessage(fr), "\n") else {
  lp <- Filter(function(l) identical(l$dpar, "mu"), fr$linpreds)[[1]]
  cat("NA in f: rows", nrow(lp$X), "spreads", env$predictor_col_spread(lp), "\n")
}
