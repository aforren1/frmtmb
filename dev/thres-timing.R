# Lane thres: what grouped thresholds cost per fit, and that a model
# without thres() is untouched.
#
#   Rscript dev/thres-timing.R lane   # the lane build
#   Rscript dev/thres-timing.R base   # the round's reference build
#
# 1. Untouched: an ordinal model without thres() fitted on each build;
#    the script prints the objective, the iteration count and the
#    estimates at full precision so the two runs can be compared
#    bitwise.
# 2. Cost (lane only): n = 6000, three groups with 4, 3 and 5
#    categories. Arms, interleaved in one process: the ungrouped model
#    y ~ x (CONTROL A and CONTROL B, the same call twice, whose ratio
#    must read 1.0), and y | thres(gr = g) ~ x. Each arm is a block of
#    fits long enough to pass 1.2 s; the minimum over 5 rounds is kept.
#    Every family is timed. Seed 21.
which <- commandArgs(trailingOnly = TRUE)[1]
lib <- if (identical(which, "base")) NULL else "/opt/rlib/lane-thres"
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(21)
n <- 6000
d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
u <- rlogis(n, 0.7 * d$x)
tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6), c = c(-1.5, -0.4, 0.5,
                                                        1.6))
d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]), 1L)

cat("\n== 1. no thres(): full-precision record ==\n")
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  ff <- get(fam, envir = asNamespace("frmtmb"))
  fit <- frm(y ~ x, data = d, family = ff())
  cat(sprintf("%-10s obj %.17g iter %d est %s\n", fam, fit$opt$objective,
              fit$opt$iterations,
              paste(sprintf("%.17g", fit$opt$par), collapse = " ")))
}
if (identical(which, "base")) quit(save = "no")

cat("\n== 2. fit time, minimum of 5 rounds, seconds per fit ==\n")
block <- function(expr_fun) {
  k <- 0L
  t0 <- proc.time()[["elapsed"]]
  repeat {
    expr_fun()
    k <- k + 1L
    el <- proc.time()[["elapsed"]] - t0
    if (el > 1.2) return(el / k)
  }
}
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  ff <- get(fam, envir = asNamespace("frmtmb"))
  arms <- list(
    controlA = function() frm(y ~ x, data = d, family = ff()),
    controlB = function() frm(y ~ x, data = d, family = ff()),
    grouped = function() frm(y | thres(gr = g) ~ x, data = d, family = ff())
  )
  res <- matrix(NA_real_, 5, length(arms), dimnames = list(NULL, names(arms)))
  for (r in 1:5) for (a in names(arms)) res[r, a] <- block(arms[[a]])
  m <- apply(res, 2, min)
  cat(sprintf(paste("%-10s controlA %.3f controlB %.3f grouped %.3f |",
                    "B/A %.2f grouped/A %.2f\n"), fam,
              m[["controlA"]], m[["controlB"]], m[["grouped"]],
              m[["controlB"]] / m[["controlA"]],
              m[["grouped"]] / m[["controlA"]]))
}
