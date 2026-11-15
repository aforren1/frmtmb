# Lane predraw2: the narrowed sides' ratio r of the fits whose draws the
# final build changed without the check firing
# (dev/predraw2-meas2-cmp.txt): each should have a side in (0.4, 0.5),
# where the ramp of draw_edge_blend() acts. Usage: Rscript <this> <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/predraw2-cases.R")
for (id in c("zip:9", "G6:1", "G6:17", "R0.1:40", "S0.1:6", "S0.1:23",
             "S0.2:22", "S0.2:23", "bglmm:26")) {
  fit <- case_fit(id)$fit
  invisible(fq(predict(fit, ndraws = 2)))
  law <- fit$cache[["draw_edge_law"]]
  vs <- c(law$lo, law$hi)
  vs <- vs[vapply(vs, length, 1L) > 1L]
  r <- vapply(vs, function(v) v[2] / 2, 0)
  cat(sprintf("%-10s r of narrowed sides: %s\n", id,
              paste(sprintf("%.3f", sort(r)), collapse = " ")))
}
