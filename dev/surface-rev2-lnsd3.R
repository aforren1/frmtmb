# Reviewer re-check, B2: map where logitnormal_sd()'s integrate() path
# errors, over p = 10^-(1..300) and logit SD 4.5 to 1e7.
source("dev/surface-rev-env.R"); rev_env("lane")
suppressPackageStartupMessages(library(frmtmb))
lnsd <- frmtmb:::logitnormal_sd
lp <- c(1, 2, 3, 5, 8, 12, 20, 30, 50, 100, 200, 300)
ss <- c(4.5, 6, 10, 20, 50, 100, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 1e6, 1e7)
bad <- 0; badmid <- 0
for (k in lp) for (s in ss) {
  p <- 10^-k
  r <- tryCatch(lnsd(p, 1 - p, s), error = function(e) NA)
  if (is.na(r)) { bad <- bad + 1; cat(sprintf("ERROR at p 1e-%d s %g\n", k, s))
    if (s <= 1000) badmid <- badmid + 1 }
}
cat("errors:", bad, "of", length(lp) * length(ss), "| with s <= 1000:", badmid, "\n")
