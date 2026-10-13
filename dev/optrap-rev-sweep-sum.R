# Reviewer of lane optrap, item 2: the lane arms rerun on the final build
# (dev/optrap-rev-log/sweep/) against the lane's own rellib-r7 arms
# (dev/optrap-log/sweep/*-r7.tsv) and the lane's earlier lane arms.
rd <- function(f) read.delim(f, quote = "", colClasses = "character")
for (v in c("spec", "mo")) for (b in c("ref", "ob26")) {
  r <- rd(sprintf("dev/optrap-log/sweep/%s-%s-r7.tsv", v, b))
  l <- rd(sprintf("dev/optrap-rev-log/sweep/%s-%s-lane.tsv", v, b))
  o <- rd(sprintf("dev/optrap-log/sweep/%s-%s-lane.tsv", v, b))
  m <- merge(r, l, by = "seed", suffixes = c(".r7", ".lane"))
  m <- merge(m, o[, c("seed", "status", "objective")], by = "seed")
  cat(sprintf(paste0("%-4s %-4s seeds %d | errors r7 %d lane(final) %d ",
                     "lane(recorded) %d | r7 errors where lane fits %d | ",
                     "lane errors where r7 fits %d | final vs recorded lane:",
                     " same status %d, same objective %d\n"),
              v, b, nrow(m), sum(m$status.r7 == "error"),
              sum(m$status.lane == "error"), sum(m$status == "error"),
              sum(m$status.r7 == "error" & m$status.lane == "ok"),
              sum(m$status.lane == "error" & m$status.r7 == "ok"),
              sum(m$status == m$status.lane),
              sum(m$objective == m$objective.lane)))
  ok <- m$status.r7 == "ok" & m$status.lane == "ok"
  dobj <- as.numeric(m$objective.lane[ok]) - as.numeric(m$objective.r7[ok])
  cat(sprintf("          both fit %d: lane lower by > 1e-6 %d, higher by > 1e-6 %d (max +%.3g)\n",
              sum(ok), sum(dobj < -1e-6, na.rm = TRUE),
              sum(dobj > 1e-6, na.rm = TRUE), max(dobj, na.rm = TRUE)))
  x <- m[m$status.r7 == "error" & m$status.lane == "ok", ]
  cat("          r7 error -> lane code:", paste(table(x$code.lane)),
      "names", paste(names(table(x$code.lane))), "| lane nwarn:",
      paste(x$nwarn.lane, collapse = ","), "\n")
}
