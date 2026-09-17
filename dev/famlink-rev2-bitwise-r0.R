## Round 1b lane against the round 0 lane build, same pairs: bitwise.
a <- readRDS("dev/famlink-rev-bitwise-lane-round0.rds"); b <- readRDS("dev/famlink-rev-bitwise-lane.rds")
k <- intersect(names(a), names(b)); nd <- 0L; nfit <- 0L
for (i in k) if (a[[i]]$status == "fit" || b[[i]]$status == "fit") { nfit <- nfit + 1L
  if (!identical(a[[i]][c("status", "ll", "est", "se", "fitted")], b[[i]][c("status", "ll", "est", "se", "fitted")])) { nd <- nd + 1L; cat(i, "differs\n") } }
cat("pairs fitted in either:", nfit, " identical round 0 vs round 1b:", nfit - nd, "\n")
