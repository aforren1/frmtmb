# Reviewer: hand-translation ok counts, reviewer build vs the lane's
# lane and r6 runs (dev/surface-bv-out/{r6,lane}).
rd <- function(dir) {
  fs <- list.files(dir, "^brms.*[.]csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
  d <- d[d$kind != "data", ]
  k <- paste(d$vignette, d$label, sep = " | ")
  d$key <- paste(k, stats::ave(seq_along(k), k, FUN = seq_along))
  d
}
x <- list(r6 = rd("dev/surface-bv-out/r6"), lane = rd("dev/surface-bv-out/lane"),
          rev = rd("dev/surface-rev-bv-out"))
cat("columns:", names(x$rev), "\n")
okcol <- if ("ok" %in% names(x$rev)) "ok" else "status"
isok <- function(d) if (okcol == "ok") as.logical(d$ok) else d$status == "ok"
for (n in names(x)) cat(n, "rows", nrow(x[[n]]), "ok", sum(isok(x[[n]])), "\n")
m <- match(x$lane$key, x$rev$key)
cat("rows matched lane->rev:", sum(!is.na(m)), "\n")
dif <- which(isok(x$lane) != isok(x$rev)[m])
cat("rows whose ok differs lane vs reviewer:", length(dif), "\n")
for (i in dif) cat("  ", x$lane$key[i], ":", isok(x$lane)[i], "->", isok(x$rev)[m[i]], "\n")
