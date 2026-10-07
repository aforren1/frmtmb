# Reviewer of lane optima, final check: compares dev/optima-rev3-nomo.R's
# four runs (base and lane, two each, pinned to core 0).
#   Rscript dev/optima-rev3-nomo-cmp.R
p <- "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev3-out/nomo-"
r <- lapply(c(b1 = "base-1", b2 = "base-2", l1 = "lane-1", l2 = "lane-2"),
            function(t) readRDS(paste0(p, t, ".rds")))
for (nm in names(r$b1)) {
  cat(nm, ": base runs identical", identical(r$b1[[nm]], r$b2[[nm]]),
      "| lane runs identical", identical(r$l1[[nm]], r$l2[[nm]]),
      "| base vs lane table identical",
      identical(r$b1[[nm]]$tab, r$l1[[nm]]$tab),
      "| messages identical", identical(r$b1[[nm]]$msg, r$l1[[nm]]$msg), "\n")
}
