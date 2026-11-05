# Compare the two outputs of review-predread-regress.R digit for digit.
a <- readRDS("dev/review-predread-out/regress-before.rds")
b <- readRDS("dev/review-predread-out/regress-after.rds")
cat("keys before", length(a), "after", length(b), "\n")
flat <- function(x) {
  if (is.character(x)) return(x)
  if (is.list(x)) return(unlist(lapply(x, flat), use.names = FALSE))
  if (is.factor(x)) return(as.character(x))
  as.numeric(as.matrix(x))
}
rows <- lapply(intersect(names(a), names(b)), function(k) {
  x <- a[[k]]
  y <- b[[k]]
  ex <- is.character(x) && length(x) == 1L && grepl("^ERROR", x)
  ey <- is.character(y) && length(y) == 1L && grepl("^ERROR", y)
  fx <- flat(x)
  fy <- flat(y)
  rel <- NA_real_
  if (is.numeric(fx) && is.numeric(fy) && length(fx) == length(fy)) {
    ok <- is.finite(fx) & is.finite(fy)
    rel <- if (any(ok)) max(abs(fx[ok] - fy[ok]) / pmax(abs(fx[ok]), abs(fy[ok]), 1e-300)) else 0
    if (any(is.na(fx) != is.na(fy))) rel <- Inf
  }
  data.frame(key = k, identical = identical(x, y) || identical(fx, fy), rel = rel,
             err_before = if (ex) substr(x, 1, 90) else "",
             err_after = if (ey) substr(y, 1, 90) else "",
             len = paste(length(fx), length(fy)),
             names_differ = !identical(names(x), names(y)))
})
r <- do.call(rbind, rows)
cat("cells", nrow(r), "bitwise", sum(r$identical), "\n")
cat("both refuse same message", sum(r$err_before != "" & r$err_before == r$err_after), "\n")
cat("rel bins: <=1e-15", sum(r$rel <= 1e-15, na.rm = TRUE), " <=1e-12", sum(r$rel <= 1e-12, na.rm = TRUE),
    " <=1e-8", sum(r$rel <= 1e-8, na.rm = TRUE), " >1e-8", sum(r$rel > 1e-8, na.rm = TRUE), "\n")
print(r[(!r$identical & (is.na(r$rel) | r$rel > 1e-12)) | r$err_before != r$err_after | r$names_differ,
        c("key", "rel", "err_before", "err_after", "len", "names_differ")], right = FALSE)
write.table(r, "dev/review-predread-out/regress-cmp.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
