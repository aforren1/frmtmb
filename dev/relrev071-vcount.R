x <- read.csv("C:/Users/adf44/source/r/frmtmb-wt-release/dev/remeasure-out/merged/port/results-merged.csv")
print(table(x$kind, x$class == "CLEAN"))
