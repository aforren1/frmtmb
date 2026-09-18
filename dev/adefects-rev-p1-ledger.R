setwd("C:/Users/adf44/source/r/frmtmb-wt-adefects")
led <- read.delim("dev/brmsport-ledger.tsv", quote = "",
                  stringsAsFactors = FALSE)
cat("ledger rows:", nrow(led), "\n")
print(table(led$outcome))
cat("\npass by class:\n"); print(table(led$class[led$outcome == "pass"]))
cat("\nany fit-data left:", sum(led$class == "fit-data"), "\n")
recs <- list.files("dev/brmsport-log", "^rec-.*[.]tsv$", full.names = TRUE)
cols <- c("kind", "pkg", "id", "verdict", "held", "vacuous", "msg",
          "caught", "raw_held")
rec <- do.call(rbind, lapply(recs, function(f) {
  x <- utils::read.delim(f, header = FALSE, quote = "", col.names = cols,
                         colClasses = "character", na.strings = NULL)
  x[x$kind == "assert", ]
}))
rec$held <- rec$held == "TRUE"
rec$raw_held <- rec$raw_held == "TRUE"
cat("
rec assert rows:", nrow(rec), "
")
w <- rec$raw_held & !rec$held
cat("withheld:", sum(w),
    " vacuous:", sum(w & startsWith(rec$msg, "VACUOUS")),
    " stale:", sum(w & startsWith(rec$msg, "STALE")),
    " parts sum to the total:",
    sum(w & startsWith(rec$msg, "VACUOUS")) +
      sum(w & startsWith(rec$msg, "STALE")) == sum(w), "
")
cat("hollow:", sum(led$class == "hollow"), "
")
