# Reviewer, claim 7: the brms_port() records of the ten core suite files,
# run gated in record mode (dev/aterms2-rev-suite.sh), compared row by
# row: lane build vs base build (both with the lane's helper), and lane
# build with the lane's helper vs with the base helper
# (dev/aterms2-rev-tt-basehelper). Log: dev/aterms2-rev-log-08-portcmp.txt
dir <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/aterms2-rev-rec"
rd <- function(tag, arm) {
  fs <- list.files(dir, sprintf("^rec-%s-%s-.*[.]tsv$", tag, arm),
                   full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    x <- utils::read.delim(f, header = FALSE, quote = "",
                           colClasses = "character")
    names(x) <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
                  "caught", "raw")[seq_len(ncol(x))]
    x$file <- basename(f)
    x
  }))
}
lane <- rd("port", "lane")
base <- rd("port", "base")
hb <- rd("porthb", "lane")
key <- function(x) paste(x$kind, x$id)
for (nm in c("lane", "base", "hb")) {
  x <- get(nm)
  cat(nm, ": rows", nrow(x), " files", length(unique(x$file)),
      " held", sum(x$held == "TRUE"), " vacuous", sum(x$vacuous == "TRUE"),
      " duplicated keys", sum(duplicated(key(x))), "\n")
}
cmp <- function(a, b, la, lb) {
  m <- merge(a[, c("kind", "id", "verdict", "held", "vacuous", "msg")],
             b[, c("kind", "id", "held", "vacuous", "msg")],
             by = c("kind", "id"), suffixes = c(".a", ".b"), all = TRUE)
  ch <- m[is.na(m$held.a) | is.na(m$held.b) | m$held.a != m$held.b |
            m$vacuous.a != m$vacuous.b, ]
  cat("\n", la, "vs", lb, ": rows compared", nrow(m), " changed", nrow(ch),
      "\n")
  for (i in seq_len(nrow(ch))) {
    cat(sprintf("  %-7s %-16s verdict=%-16s held %s->%s vacuous %s->%s\n",
                ch$kind[i], ch$id[i], ch$verdict[i], ch$held.b[i],
                ch$held.a[i], ch$vacuous.b[i], ch$vacuous.a[i]))
    cat("      ", lb, ":", substr(ch$msg.b[i], 1, 150), "\n")
    cat("      ", la, ":", substr(ch$msg.a[i], 1, 150), "\n")
  }
}
cmp(lane, base, "lane", "base")
cmp(lane, hb, "lane-helper", "base-helper")
