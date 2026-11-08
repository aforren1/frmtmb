# Lane arstep: what the second descent did, from the traced runs
# (dev/arstep-ng3.R, dev/arstep-r2ng3.R; traces in dev/arstep-trace.R).
# Per set: rows, rows with a second descent, rows it raised by more than
# 1e-3, the seconds it took against the fits' total, and the working
# marginal's gap from its best to the second maximum on rows that gained
# and rows that did not. Args: none.
out <- "C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-out"
one <- function(f, tag, sdtag) {
  l <- readLines(file.path(out, f))
  fits <- grep(paste0("^", tag, " "), l, value = TRUE)
  sd <- grep(paste0("^", sdtag, " "), l, value = TRUE)
  wall <- if (tag == "ANG") as.numeric(sub(".* wall ([^ ]+).*", "\\1", fits)) else
    as.numeric(sub(".* ([0-9.]+)s$", "\\1", fits))
  ran <- grepl(" ran first ", sd)
  num <- function(s, key) as.numeric(sub(paste0(".* ", key, " ([^ ;]+).*"), "\\1", s))
  first <- num(sd[ran], "first"); second <- num(sd[ran], "second"); secs <- num(sd[ran], "secs")
  gain <- second - first
  mb <- num(sd[ran], "best")
  ma <- as.numeric(sub(",.*", "", sub(".* alts ([^;]+);.*", "\\1", sd[ran])))
  dm <- ma - mb
  key <- sub(paste0("^", sdtag, " [^ ]+ "), "", sub(" marginal.*| none$| ran .*", "", sd))
  cat(sprintf("TSUM %-6s rows %d, second descent on %d, raised by >1e-3 on %d, lowered (kept the first) on %d; seconds %.1f of %.1f total\n",
              f, length(fits), sum(ran), sum(gain > 1e-3), sum(gain < -1e-3), sum(secs), sum(wall)))
  cat(sprintf("TSUM %-6s marginal gap to the second maximum: rows that gained %s; others: median %.2f, range %.2f to %.2f\n", f,
              paste(sprintf("%s %.2f (+%.3g)", key[ran][gain > 1e-3], dm[gain > 1e-3], gain[gain > 1e-3]), collapse = "; "),
              median(dm[gain <= 1e-3]), min(dm[gain <= 1e-3]), max(dm[gain <= 1e-3])))
}
one("ng3-lane.txt", "ANG", "ANGSD")
one("r2ng3-lane.txt", "RNG2", "RNG2SD")
cat("DONE\n")
