# The shipped documentation's numbers, against the lane's own record.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
cat("\n=== the vignette's '0.2 percent short' ===\n")
P <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1, V = 10)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                 addl = 0L, ss = TRUE)
ts <- c(0, 1, 4, 8)
a <- frm_lincmt(parms = P, times = ts, ncmt = 2, depot = TRUE,
                events = ev, n_ss = 20L)
b <- frm_lincmt(parms = P, times = ts, ncmt = 2, depot = TRUE,
                events = ev)
cat("  n_ss = 20 against the exact limit, on the schedule the",
    "vignette names:\n")
cat("   max |a - b| / max|b| =", format(max(abs(a - b)) / max(abs(b)),
                                        digits = 4), "\n")
cat("   findings.md says 3.96e-03; the vignette says 0.2 percent\n")
