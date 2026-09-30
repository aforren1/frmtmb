# Reviewer helper (lane ceplot): every mutant pattern matches its
# installed function exactly once.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-mutants.R")
for (nm in names(mutants)) {
  for (ed in mutants[[nm]]) {
    txt <- paste(deparse(get(ed[[2]], asNamespace(ed[[1]]))), collapse = "\n")
    n <- sum(gregexpr(ed[[3]], txt, fixed = TRUE)[[1]] > 0)
    cat(sprintf("%-18s %-32s matches %d\n", nm, ed[[2]], n))
  }
}
