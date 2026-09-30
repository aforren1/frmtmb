# Sourced by postfit2-rev-p2-mmsplit.R: give each new member its own
# draw-cache key (deparse splits this call across lines).
tx <- paste(tx, collapse = "\n")
tx2 <- sub("paste0[(]\"new:\",[[:space:]]*u[[]j[]][)]",
           "paste0(\"new:\", u[j], \"#\", j)", tx)
stopifnot(!identical(tx, tx2))
tx <- tx2
