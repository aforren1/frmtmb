# Reviewer, lane wt-conditions (recheck 1): compare the base and lane
# records of dev/conditions-rev-falsealarm.R case by case.
#   Rscript dev/conditions-rev-falsealarm-cmp.R [suffix]
sfx <- commandArgs(TRUE)[1L]; if (is.na(sfx)) sfx <- ""
o <- "C:/Users/adf44/source/r/frmtmb-wt-conditions/dev/conditions-rev-log"
b <- readRDS(file.path(o, paste0("falsealarm", sfx, "-base.rds")))$cases
l <- readRDS(file.path(o, paste0("falsealarm", sfx, "-lane.rds")))$cases
m <- merge(b, l, by = "case", suffixes = c(".base", ".lane"))
ok_b <- m$result.base == "OK"; ok_l <- m$result.lane == "OK"
cat("cases", nrow(m), " OK both", sum(ok_b & ok_l), " error both",
    sum(!ok_b & !ok_l), " REGRESSION (base OK, lane error)",
    sum(ok_b & !ok_l), " lane OK, base error", sum(!ok_b & ok_l), "\n\n")
show <- function(sel, title) {
  cat("==", title, "\n")
  for (i in which(sel)) {
    cat(sprintf("%s\n   base: %s\n   lane: %s\n", m$case[i],
                substr(m$result.base[i], 1, 200),
                substr(m$result.lane[i], 1, 200)))
  }
  cat("\n")
}
show(ok_b & !ok_l, "REGRESSIONS")
show(!ok_b & ok_l, "lane OK, base error")
show(!ok_b & !ok_l, "error on both")
