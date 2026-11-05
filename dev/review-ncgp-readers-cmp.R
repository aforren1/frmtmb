# Compares the two RDS files of review-ncgp-readers.R: rows of largest relative
# difference, every refusal or NaN that exists on one build only.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
a <- readRDS(file.path(REV_OUT, "readers-lane.rds"))
b <- readRDS(file.path(REV_OUT, "readers-base.rds"))
fl <- function(v) { if (is.data.frame(v)) v <- as.matrix(v); suppressWarnings(as.numeric(unlist(v))) }
rows <- list()
for (nm in names(a)) {
  if (!is.list(a[[nm]]) || !is.list(b[[nm]]) || is.null(a[[nm]]$ll) || is.null(b[[nm]]$ll)) { cat("RR", nm, "FIT lane:", a[[nm]]$fit %||% "ok", "base:", b[[nm]]$fit %||% "ok", "\n"); next }
  for (k in names(a[[nm]])) {
    x <- a[[nm]][[k]]; y <- b[[nm]][[k]]
    ex <- is.character(x) && length(x) == 1 && grepl("^ERR", x); ey <- is.character(y) && length(y) == 1 && grepl("^ERR", y)
    if (ex || ey) {
      if (ex && ey) { rows[[length(rows) + 1L]] <- data.frame(m = nm, k = k, s = if (identical(x, y)) "both-refuse" else "DIFF-REFUSE", d = NA, x = substr(x, 6, 70)); next }
      rows[[length(rows) + 1L]] <- data.frame(m = nm, k = k, s = if (ex) "LANE-REFUSES" else "BASE-REFUSES", d = NA, x = substr(if (ex) x else y, 6, 90)); next
    }
    fx <- fl(x); fy <- fl(y)
    if (length(fx) != length(fy)) { rows[[length(rows) + 1L]] <- data.frame(m = nm, k = k, s = sprintf("LEN %d %d", length(fx), length(fy)), d = NA, x = ""); next }
    na <- is.na(fx) != is.na(fy)
    ok <- !is.na(fx) & !is.na(fy)
    d <- if (any(ok)) max(abs(fx[ok] - fy[ok]) / pmax(abs(fy[ok]), 1e-3)) else NA
    rows[[length(rows) + 1L]] <- data.frame(m = nm, k = k, s = if (any(na)) sprintf("NA %d lane %d base", sum(is.na(fx)), sum(is.na(fy))) else "ok", d = d, x = "")
  }
}
r <- do.call(rbind, rows)
cat("rows", nrow(r), "\n")
print(table(r$s))
print(r[r$s != "ok" & r$s != "both-refuse", ], row.names = FALSE)
print(utils::head(r[order(-r$d), c("m", "k", "d")], 45), row.names = FALSE)
print(tapply(r$d, r$k, max, na.rm = TRUE))
