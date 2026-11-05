# Short check of fix round 1, item 3: the message's cost at 450 positions. Seed 6602.
#   Rscript dev/review-gpgap-cost.R
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(6602)
n <- 450
d <- data.frame(x = round(runif(n, 0, 6), 4))
d$y <- sin(2 * d$x) + 0.1 * rnorm(n)
t1 <- proc.time()[["elapsed"]]
f <- suppressWarnings(frm(bf(y ~ gp(x)), data = d))
t2 <- proc.time()[["elapsed"]]
m <- frmtmb:::gp_nugget_message(f)
t3 <- proc.time()[["elapsed"]]
cat(sprintf("COST n %d fit %.1f s, message %.2f s (%.2f percent), fires %s\n", n, t2 - t1,
            t3 - t2, 100 * (t3 - t2) / (t2 - t1), !is.null(m)))
