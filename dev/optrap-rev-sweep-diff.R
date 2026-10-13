# Reviewer of lane optrap, item 2: seeds where both builds fit and the
# objectives differ by more than 1e-6.
rd <- function(f) read.delim(f, quote = "", colClasses = "character")
for (v in c("spec", "mo")) for (b in c("ref", "ob26")) {
  r <- rd(sprintf("dev/optrap-log/sweep/%s-%s-r7.tsv", v, b))
  l <- rd(sprintf("dev/optrap-rev-log/sweep/%s-%s-lane.tsv", v, b))
  o <- rd(sprintf("dev/optrap-log/sweep/%s-%s-lane.tsv", v, b))
  m <- merge(r, l, by = "seed", suffixes = c(".r7", ".lane"))
  m$rec <- o$objective[match(m$seed, o$seed)]
  ok <- m$status.r7 == "ok" & m$status.lane == "ok"
  m$d <- as.numeric(m$objective.lane) - as.numeric(m$objective.r7)
  x <- m[ok & abs(m$d) > 1e-6, c("seed", "code.r7", "objective.r7",
                                 "code.lane", "objective.lane", "rec", "d")]
  cat("==", v, b, "| final lane objective == recorded lane objective:",
      sum(m$objective.lane == m$rec, na.rm = TRUE), "of",
      sum(!is.na(m$rec) & m$rec != "NA"), "\n")
  print(x)
}
