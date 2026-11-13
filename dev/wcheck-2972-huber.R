# Lane wcheck, item 2972: huber() reached a log likelihood of +1.39e5 on
# 400 rows on rellib-r8 (dev/optstart-battery.R, design hir2, seeds 2,
# 3, 4, 7). This is that design for huber alone, copied from the
# battery: 400 rows, 20 groups, slope 1e3 times sigma = 1, group sd 0.5,
# y ~ x + (1 | g), unweighted. Each fit from the default start and from
# the true values.
#   Rscript dev/wcheck-2972-huber.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
one <- function(d, start = NULL) {
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(
    frm(y ~ x + (1 | g), data = d, family = huber(), start = start)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(f, "error")) return(sprintf("ERROR %s", conditionMessage(f)))
  sprintf("logLik %.6g code %d warnings %d", as.numeric(logLik(f)),
          f$opt$convergence, length(w))
}
for (s in 1:10) {
  set.seed(s)
  n <- 400; G <- 20
  g <- factor(rep(seq_len(G), length.out = n))
  x <- rnorm(n)
  b0 <- 0.3; b1 <- 1e3; sdg <- 0.5
  bg <- rnorm(G, 0, sdg)
  mu <- b0 + b1 * x + bg[g]
  y <- mu + rnorm(n, 0, 1)
  d <- data.frame(y = y, x = x, g = g)
  truth <- list(beta = c(`(Intercept)` = b0, x = b1),
                betad = c(`sigma_(Intercept)` = 0), theta = log(sdg))
  cat(sprintf("seed %2d default: %s | truth start: %s\n", s, one(d),
              one(d, truth)))
}
cat("DONE\n")
