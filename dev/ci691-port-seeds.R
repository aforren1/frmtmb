# Lane ci691: the pin of test-optima.R ("a finite trial PORT rejected is
# not reported either") in a fresh process: run_optimizer() on the cs()
# ordinal mixture, seeds 1, 9 and 13, default start, and whether
# fn(par) is the reported objective. Run in several processes to see
# that the outcome does not follow the core type (rtmb-pitfalls item 21).
#   Rscript dev/ci691-port-seeds.R <lib or base>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r7"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressMessages(library(frmtmb))
optima_csmix_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
out <- character(0)
for (s in c(1, 9, 13)) {
  u <- suppressWarnings(frm(bf(y ~ x + cs(z)),
                            family = mixture(cumulative(), sratio()),
                            data = optima_csmix_data(s),
                            dry_run = "objective"))
  obj <- u$obj
  r <- suppressWarnings(frmtmb:::run_optimizer(
    "nlminb", obj$par, obj$fn, obj$gr, -Inf, Inf,
    list(eval.max = 1000, iter.max = 1000)))
  out <- c(out, sprintf("seed %d %s", s,
                        if (identical(obj$fn(r$par), r$objective)) "same"
                        else "DIFF"))
}
cat("lib", find.package("frmtmb"), ":", paste(out, collapse = ", "), "\n")
