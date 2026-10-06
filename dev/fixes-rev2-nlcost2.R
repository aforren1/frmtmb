# Reviewer of lane fixes, re-check: how check_nl_identified()'s time
# grows with n, against the fit of the same model (one fit per n), and
# the guard against a no-op control on the same frame.
#   Rscript dev/fixes-rev2-nlcost2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
guard <- get("check_nl_identified", asNamespace("frmtmb"))
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
fo <- bf(y ~ a * exp(b * x) + c0, a ~ 1 + z, b ~ 1 + z, c0 ~ 1, nl = TRUE)
for (n in c(2000, 8000, 32000)) {
  set.seed(31)
  d <- data.frame(x = runif(n), z = rnorm(n))
  d$y <- 2 * exp((0.5 + 0.2 * d$z) * d$x) + 1 + rnorm(n, 0, 0.3)
  fr <- suppressMessages(frm(fo, data = d, dry_run = "frame"))
  g <- tm(guard(fr$spec, fr, NULL))
  ctl <- tm(invisible(NULL))
  ft <- tm(suppressMessages(suppressWarnings(
    frm(fo, data = d, start = list(beta = c(2, 0, 0.5, 0, 1))))))
  cat(sprintf("n = %6d: guard %.2f s, control %.2f s, fit (with guard) %.2f s\n",
              n, g, ctl, ft))
}
