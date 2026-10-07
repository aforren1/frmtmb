# Compare fragility-scan runs (dev/ciharden-scan.sh) expectation by
# expectation and list what depends on the configuration.
#
#   Rscript dev/ciharden-scan-diff.R <out prefix> <run> <run> [<run> ...]
#
# A "signature" is one non-passing expectation: kind (FAIL, ERROR,
# WARN, SKIP), file:line, test name and message with every number
# masked, so a failure whose printed numbers move is still one
# signature. A signature present in some runs and absent in others is
# configuration-dependent; one present in every run is a standing
# failure or skip, listed apart. A test whose pass count differs
# between runs is listed too, since a loop over what a fit returns can
# move it without any failure. Files without a RESULT line are listed
# first, because their counts mean nothing.
#
# Writes <out prefix>.tsv (one row per varying signature or count, with
# a 0/1 or count column per run) and prints a summary. Before calling a
# row a BLAS effect, compare it with two runs of one configuration:
# fresh processes on this hybrid CPU can differ by a few ulps
# (dev/rtmb-pitfalls.md item 21), and wall-clock tests move with load.
a <- commandArgs(TRUE)
out <- a[1]
runs <- a[-1]
logd <- file.path(dirname(dirname(normalizePath(
  sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))),
  "dev", "ciharden-log")
mask <- function(s) {
  s <- gsub("0x[0-9a-fA-F]+", "#", s)
  s <- gsub("-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?", "#", s)
  gsub("[ ]+", " ", s)
}
read_run <- function(r) {
  d <- file.path(logd, paste0(r, "-suite"))
  fs <- list.files(d, pattern = "[.]txt$", full.names = TRUE)
  sig <- list()
  cnt <- list()
  res <- list()
  for (f in fs) {
    key <- sub("[.]txt$", "", basename(f))
    l <- readLines(f, warn = FALSE)
    res[[key]] <- any(grepl("^RESULT .* pass=", l))
    det <- grep("^DETAIL\t", l, value = TRUE)
    if (length(det)) {
      p <- strsplit(det, "\t", fixed = TRUE)
      s <- vapply(p, function(x) {
        paste(key, x[2], x[3], x[4], mask(x[5]), sep = "\t")
      }, "")
      sig[[key]] <- unique(s)
    }
    tl <- grep("^TEST\t", l, value = TRUE)
    if (length(tl)) {
      p <- strsplit(tl, "\t", fixed = TRUE)
      cnt[[key]] <- stats::setNames(
        as.integer(vapply(p, `[`, "", 3)),
        paste(key, vapply(p, `[`, "", 2), sep = "\t"))
    }
  }
  list(sig = unlist(sig, use.names = FALSE), cnt = unlist(unname(cnt)),
       res = unlist(res))
}
R <- lapply(runs, read_run)
names(R) <- runs
files <- sort(unique(unlist(lapply(R, function(x) names(x$res)))))
cat("runs:", paste(runs, collapse = " "), "\n")
for (r in runs) {
  cat(sprintf("%s: %d files, %d without a RESULT line\n", r,
              length(R[[r]]$res), sum(!R[[r]]$res)))
  bad <- names(R[[r]]$res)[!R[[r]]$res]
  if (length(bad)) cat("  no RESULT:", paste(bad, collapse = " "), "\n")
}
allsig <- sort(unique(as.character(unlist(lapply(R, `[[`, "sig")))))
M <- vapply(R, function(x) as.integer(allsig %in% x$sig),
            integer(length(allsig)))
if (!is.matrix(M)) M <- matrix(M, ncol = length(runs))
colnames(M) <- runs
vary <- rowSums(M) > 0 & rowSums(M) < length(runs)
always <- rowSums(M) == length(runs)
split_sig <- function(s) {
  if (!length(s)) return(NULL)
  do.call(rbind, lapply(strsplit(s, "\t", fixed = TRUE), function(x) {
    data.frame(file = x[1], kind = x[2], line = x[3], test = x[4],
               msg = substr(x[5], 1, 160))
  }))
}
rows <- list()
if (any(vary)) {
  rows[[1]] <- cbind(split_sig(allsig[vary]), what = "signature",
                     as.data.frame(M[vary, , drop = FALSE]))
}
allt <- sort(unique(unlist(lapply(R, function(x) names(x$cnt)))))
C <- vapply(R, function(x) {
  v <- x$cnt[allt]
  ifelse(is.na(v), -1L, v)
}, integer(length(allt)))
if (!is.matrix(C)) C <- matrix(C, ncol = length(runs))
colnames(C) <- runs
cv <- apply(C, 1, function(v) length(unique(v)) > 1)
if (any(cv)) {
  p <- strsplit(allt[cv], "\t", fixed = TRUE)
  rows[[2]] <- cbind(data.frame(file = vapply(p, `[`, "", 1), kind = "PASSES",
                                line = "", test = vapply(p, `[`, "", 2),
                                msg = "pass count differs (-1: not run)"),
                     what = "count", as.data.frame(C[cv, , drop = FALSE]))
}
tab <- do.call(rbind, rows)
if (is.null(tab)) {
  cat("no configuration-dependent signature or count\n")
  tab <- data.frame()
} else {
  tab <- tab[order(tab$file, tab$test, tab$kind), ]
  cat(sprintf("configuration-dependent: %d signatures, %d pass counts\n",
              sum(tab$what == "signature"), sum(tab$what == "count")))
  print(tab[, c("file", "kind", "line", runs)], row.names = FALSE)
}
utils::write.table(tab, paste0(out, ".tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)
st <- split_sig(allsig[always])
if (!is.null(st)) {
  st <- st[st$kind != "SKIP", , drop = FALSE]
  cat(sprintf("in every run, not a skip: %d\n", nrow(st)))
  if (nrow(st)) print(st[, c("file", "kind", "line", "msg")],
                      row.names = FALSE)
}
if (is.null(st)) cat("in every run, not a skip: 0\n")
