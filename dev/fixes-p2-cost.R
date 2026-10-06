# Lane fixes, punch round 2: the cost of the fitted-point flat check.
# The reviewer's model (dev/fixes-rev2-nlcost2.R), a * exp(b * x) + c0
# with a, b ~ 1 + z: per n, the whole frm() with the check and with it
# replaced by a no-op (the control), interleaved in one process, and
# the check alone on the fitted objective; the minimum over rounds.
#   Rscript dev/fixes-p2-cost.R <lib> <rounds>
a <- commandArgs(TRUE)
lib <- a[1]
rounds <- if (length(a) > 1L) as.integer(a[2]) else 3L
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
real <- get("nl_flat_message", ns)
noop <- function(obj, opt, frame) NULL
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
fo <- bf(y ~ a * exp(b * x) + c0, a ~ 1 + z, b ~ 1 + z, c0 ~ 1, nl = TRUE)
fit_it <- function(d) {
  suppressMessages(suppressWarnings(
    frm(fo, data = d, start = list(beta = c(2, 0, 0.5, 0, 1)))))
}
for (n in c(2000, 8000, 32000, 200000)) {
  set.seed(31)
  d <- data.frame(x = runif(n), z = rnorm(n))
  d$y <- 2 * exp((0.5 + 0.2 * d$z) * d$x) + 1 + rnorm(n, 0, 0.3)
  with_check <- without <- alone <- numeric(0)
  for (r in seq_len(rounds)) {
    assignInNamespace("nl_flat_message", real, ns = "frmtmb")
    with_check <- c(with_check, tm(f <- fit_it(d)))
    assignInNamespace("nl_flat_message", noop, ns = "frmtmb")
    without <- c(without, tm(fit_it(d)))
    alone <- c(alone, tm(real(f$obj, f$opt, f$frame)))
  }
  assignInNamespace("nl_flat_message", real, ns = "frmtmb")
  cat(sprintf(paste0("n = %6d: fit with check %.2f s, without (control) ",
                     "%.2f s, check alone %.3f s (min of %d rounds)\n"),
              n, min(with_check), min(without), min(alone), rounds))
}
