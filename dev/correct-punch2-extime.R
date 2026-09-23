# Does this lane make the `?residuals.frmtmb_fit` example slower, or is
# the `R CMD check` examples NOTE this box's load?
#
# The NOTE cannot answer it. The same example measured 3.05, 5.32, 10.28
# and 11.39 seconds of USER time over four passes of IDENTICAL code in
# two checks of the same tree, which is the swing dev/lane-rules.md
# records for this machine. So the arms have to be interleaved and
# carried against a control that must report 1.0.
#
# One arm per PROCESS, because two builds of one package cannot live in
# one R session. Rounds alternate base, lane, base, lane, so a load
# excursion hits both arms. Each process also runs a fixed arithmetic
# CONTROL that touches no frmtmb code.
#
# READ THE LEVEL, NOT THE RATIO. Alternation removes staleness and
# ordering bias; it does not remove variance. Six rounds here gave an
# example ratio of 0.816 with the control at 0.992, and nine rounds in
# the recheck gave 1.072 with the control at 1.048: opposite directions,
# so this does not resolve a difference between the arms at these
# counts. What it does settle is that BOTH arms run far over the 5
# second threshold, which is the only thing the NOTE turns on.
#
#   Rscript dev/correct-punch2-extime.R <base|lane> <rounds>
#   Rscript dev/correct-punch2-extime.R report
arm <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-log"
TSV <- file.path(OUT, "punch2-extime.tsv")

if (identical(arm, "report")) {
  x <- utils::read.delim(TSV)
  agg <- function(a, w) min(x$secs[x$arm == a & x$what == w])
  for (w in c("example", "control")) {
    b <- agg("base", w)
    l <- agg("lane", w)
    cat(sprintf("%-8s base %.3f  lane %.3f  lane/base %.3f\n",
                w, b, l, l / b))
  }
  cat("rounds per arm:", sum(x$arm == "base" & x$what == "example"), "\n")
  cat("all example times, in the order they ran:\n")
  e <- x[x$what == "example", ]
  cat(sprintf("  %-5s %s\n", e$arm, sprintf("%.3f", e$secs)), sep = "")
  quit(save = "no")
}

source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
if (identical(arm, "base")) .libPaths(.libPaths()[-1L])
suppressMessages(library(frmtmb))
rounds <- as.integer(commandArgs(trailingOnly = TRUE)[2])

# the example verbatim from man/residuals.frmtmb_fit.Rd, minus the plot
example_block <- function() {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  head(residuals(fit))
  head(residuals(fit, type = "pearson"))
  pr <- residuals(fit, type = "pearson")[, "Estimate"]
  sum(pr^2) / df.residual(fit)
  r <- residuals(fit, type = "osa")[, "Estimate"]
  invisible(r)
}

# a control that reaches no frmtmb code, sized past the 1.2 s the timing
# rule asks for so the 10 ms clock tick is not the measurement
control_block <- function() {
  s <- 0
  for (i in 1:220) {
    m <- matrix(stats::rnorm(200 * 200), 200)
    s <- s + sum(diag(crossprod(m)))
  }
  invisible(s)
}

rows <- list()
for (k in seq_len(rounds)) {
  for (w in c("example", "control")) {
    t0 <- proc.time()[["user.self"]]
    if (identical(w, "example")) example_block() else control_block()
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm, round = k, what = w,
      secs = proc.time()[["user.self"]] - t0)
  }
}
out <- do.call(rbind, rows)
utils::write.table(out, TSV, sep = "\t", row.names = FALSE,
                   col.names = !file.exists(TSV), append = file.exists(TSV),
                   quote = FALSE)
cat("arm", arm, "frmtmb from", dirname(find.package("frmtmb")), "\n")
print(out)
