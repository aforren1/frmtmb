# Offline comparison of probe strategies from dev/nlstart-patterns.R output,
# against the best logLik any arm of the scan reached on each row.
# Usage: Rscript dev/nlstart-patsum.R <patdir> <scandir> [<scandir> ...]
a <- commandArgs(trailingOnly = TRUE)
l <- unlist(lapply(list.files(a[1], "[.]txt$", full.names = TRUE),
                   function(f) grep("^PAT ", readLines(f), value = TRUE)))
sc <- unlist(lapply(a[-1], function(d) {
  unlist(lapply(list.files(d, "[.]txt$", full.names = TRUE),
                function(f) grep("^ROW ", readLines(f), value = TRUE)))
}))
s <- read.table(text = sub("^ROW ", "", sc), stringsAsFactors = FALSE)
skey <- paste(s$V1, s$V2, s$V3, s$V4)
best_scan <- tapply(pmax(s$V5, s$V6, na.rm = TRUE), skey, max, na.rm = TRUE)
rows <- lapply(strsplit(sub("^PAT ", "", l), " "), function(v) {
  pats <- do.call(rbind, strsplit(v[-(1:7)], ":"))
  list(key = paste(v[1], v[2], v[3], v[4]), fit = as.numeric(v[5]),
       t_fit = as.numeric(v[6]), sign = v[7], pat = pats[, 1],
       ll = as.numeric(pats[, 2]), t = as.numeric(pats[, 3]))
})
agree <- function(u, v) !is.na(u) && abs(u - v) <= 1e-6 * max(1, abs(v))
strategies <- list(
  none = function(r) integer(0),
  all = function(r) seq_along(r$pat),
  neg_first = function(r) {
    k <- length(r$pat)
    neg <- k
    if (agree(r$ll[neg], r$fit)) neg else seq_len(k)
  },
  pos_neg_first = function(r) {
    k <- length(r$pat)
    two <- unique(c(1L, k))
    if (all(vapply(two, function(i) agree(r$ll[i], r$fit), NA))) two else
      seq_len(k)
  },
  other_orthant = function(r) which(r$pat != r$sign))
for (nm in names(strategies)) {
  short <- 0L
  silent_rows <- character(0)
  t_all <- 0
  t_fit <- 0
  for (r in rows) {
    use <- strategies[[nm]](r)
    got <- max(c(r$fit, r$ll[use]), na.rm = TRUE)
    b <- max(best_scan[[r$key]], got, r$ll, na.rm = TRUE)
    if (b - got > 1e-4) {
      short <- short + 1L
      silent_rows <- c(silent_rows, r$key)
    }
    t_all <- t_all + r$t_fit + sum(r$t[use])
    t_fit <- t_fit + r$t_fit
  }
  cat(sprintf("%-14s short %3d  time ratio %.3f  %s\n", nm, short,
              t_all / t_fit, paste(head(silent_rows, 6), collapse = "; ")))
}
cat("patterns per row:", paste(names(table(sapply(rows, function(r) {
  length(r$pat)
}))), table(sapply(rows, function(r) length(r$pat))), sep = "x",
collapse = " "), "\n")
