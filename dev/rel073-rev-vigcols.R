# Release reviewer, 0.73.0: the vignette re-measure's outcome columns,
# r073 against 0.72.0's recorded r072, row by row.
#   Rscript dev/rel073-rev-vigcols.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-release")
h <- read.delim("dev/brms-port/outcomes-hand.tsv", colClasses = "character", quote = "")
m <- read.delim("dev/brms-port/outcomes-mech.tsv", colClasses = "character", quote = "")
cat("columns:", grep("r07", names(h), value = TRUE), "|",
    grep("r07", names(m), value = TRUE), "\n")
cat("hand r073.ok == r072.ok:", sum(h$r073.ok == h$r072.ok), "of", nrow(h), "\n")
for (k in c("raw", "spell", "class")) {
  cat("mech", k, sum(m[[paste0("r073.", k)]] == m[[paste0("r072.", k)]]),
      "of", nrow(m), "\n")
}
