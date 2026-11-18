# Reviewer short check: cum_sratio n 1000 seed 2, w = 1/1000 ends 0.226
# (per unit) below w = 1 after the swap. Is it the refit from the swapped
# start (the optimizer under a scaled objective) or the settle? nlminb
# from the same swapped start on each weight's objective, and the
# default's record.
#   Rscript ordstarts-short-p10.R > ordstarts-short-p10.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(2)
n <- 1000
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
q <- function(e) suppressWarnings(suppressMessages(e))
fam <- mixture(cumulative(), sratio())
for (w in c(1, 1 / 1000)) {
  d$w <- w
  a0 <- q(frm(bf(y | weights(w) ~ x + z), family = fam, data = d,
              control = frmtmb_control(starts = 0)))
  D <- q(frm(bf(y | weights(w) ~ x + z), family = fam, data = d))
  sw <- ns$mixture_ord_swap_starts(a0, "y")[[1]]
  o <- nlminb(sw$par, a0$obj$fn, a0$obj$gr,
              control = list(eval.max = 5000, iter.max = 5000))
  ms <- D$opt$mixture_starts
  cat(sprintf(paste0("w %-6g | A0 ll/w %.4f | D ll/w %.4f (code %d, %s) |",
                     " nlminb from the swapped start: -obj/w %.4f (conv %d: %s)\n"),
              w, as.numeric(logLik(a0)) / w, as.numeric(logLik(D)) / w,
              D$opt$convergence,
              paste(names(ms), signif(ms, 4), sep = "=", collapse = " "),
              -o$objective / w, o$convergence, o$message))
}
