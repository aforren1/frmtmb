# What frmtmb 0.69.0 (rellib-r7) says on mo(x) * f, mo(x):f and
# (mo(x) | g). Usage: Rscript dev/moparity-refusals.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
set.seed(1)
n <- 120
dd <- data.frame(
  x = sample(0:3, n, TRUE),
  f = factor(sample(c("a", "b", "c"), n, TRUE)),
  g = factor(sample(paste0("g", 1:8), n, TRUE)),
  z = rnorm(n)
)
dd$y <- rnorm(n)
try_fit <- function(fo) {
  cat("\n==", deparse1(fo), "\n")
  r <- tryCatch(frm(bf(fo), family = gaussian(), data = dd),
                error = function(e) e)
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(round(r$opt$par, 4))
}
try_fit(y ~ mo(x) * f)
try_fit(y ~ mo(x):f)
try_fit(y ~ mo(x) + (mo(x) | g))
try_fit(y ~ (mo(x) | g))
