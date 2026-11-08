# Wall time of the sign probe on the 450 nl rows of dev/nlstart-scan.R,
# measured in one process: per row the default fit with the probe off,
# on, and off again, in a rotating order. off2 / off is the control and
# must read 1. Also records each arm's logLik.
# Usage: NLSTART_ARM=after Rscript dev/nlstart-cost.R <family>
fam <- commandArgs(trailingOnly = TRUE)[1]
source("dev/nlstart-common.R")
src <- readLines("dev/nlstart-scan.R")
eval(parse(text = src[grep("^llv <-", src):(grep("^cat[(]\"ROWHEAD", src) - 1L)]))
arm_fit <- function(fm, d, probe) {
  old <- options(frmtmb.nl_sign_probe = probe)
  on.exit(options(old))
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(fm, family = F[[fam]], data = d))
  c(ll = llv(r), t = proc.time()[["elapsed"]] - t0)
}
i <- 0L
for (seed in 1:5) for (slope in c(0.6, -0.6)) for (bn in names(bodies)) {
  b <- bodies[[bn]]
  d <- gen(seed * 100 + 1, slope)
  fm <- do.call(bf, c(list(as.formula(paste("y ~", deparse1(b$f)))),
                      lapply(b$p, function(p) as.formula(paste(p, "~ 1"))),
                      list(nl = TRUE)))
  arms <- list(off = FALSE, on = TRUE, off2 = FALSE)
  ord <- c("off", "on", "off2")
  ord <- ord[(seq_along(ord) + i) %% 3L + 1L]
  i <- i + 1L
  res <- list()
  for (a in ord) res[[a]] <- arm_fit(fm, d, arms[[a]])
  cat(sprintf("COST %s %d %.1f %s %.6f %.6f %.6f %.3f %.3f %.3f\n", fam, seed,
              slope, bn, res$off[["ll"]], res$on[["ll"]], res$off2[["ll"]],
              res$off[["t"]], res$on[["t"]], res$off2[["t"]]))
  flush.console()
}
