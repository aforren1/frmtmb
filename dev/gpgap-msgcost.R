# What the nugget warning's reading costs a fit: gp_nugget_message() against
# the fit itself, 200 positions, interleaved, minimum of 5 rounds. The control is
# a second timing of the same call and must read about 1. Seed 94.
#   Rscript dev/gpgap-msgcost.R
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(94)
x <- round(rnorm(200), 3)
d <- data.frame(x = x, y = sin(2 * x) + rnorm(200, 0, 0.3))
fit <- NULL
ft <- numeric(5); m1 <- numeric(5); m2 <- numeric(5)
tm <- function(expr) {
  # grown past 1.2 s, since proc.time() ticks at 10 ms here
  k <- 1L
  repeat {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(k)) eval.parent(substitute(expr))
    t1 <- proc.time()[["elapsed"]] - t0
    if (t1 > 1.2) return(t1 / k)
    k <- k * 2L
  }
}
for (r in 1:5) {
  ft[r] <- tm(fit <- frm(bf(y ~ gp(x)), data = d))
  m1[r] <- tm(frmtmb:::gp_nugget_message(fit))
  m2[r] <- tm(frmtmb:::gp_nugget_message(fit))
}
cat(sprintf("GPCOST fit %.3f s, message %.4f s (%.2f percent of the fit), control %.2f\n",
            min(ft), min(m1), 100 * min(m1) / min(ft), min(m1) / min(m2)))
