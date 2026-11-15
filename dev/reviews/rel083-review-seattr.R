# Release review of 0.83.0: attribute the SE records that moved between
# rellib-r21 and rellib-r22 (dev/rel083-log/se-cmp.txt) to lanes. Per
# file and test whose records differ, the record sequence of each lane
# library's trace (dev/reviews/rel083-review-serun.sh) is compared with
# the merged build's and with 0.82.0's (equal: every record the same
# length and NaN pattern and within 1e-8 relative).
#   Rscript dev/reviews/rel083-review-seattr.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-rel083/dev"
lanes <- c("selimit2", "quadseam2", "slopescale", "nlstop2", "priors2",
           "famend", "silent2", "arspat2", "refits", "parity2", "wcheck2",
           "gpfield", "predraw2", "hygiene", "ci082")
num <- function(s) suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1]]))
bytest <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- readLines(f, warn = FALSE)
  test <- "(file level)"; out <- list()
  for (l in x) {
    if (startsWith(l, "TEST ")) { test <- sub("^TEST ", "", l); next }
    tag <- sub(" .*", "", l)
    if (!tag %in% c("COV", "VCT")) next
    v <- if (tag == "COV") sub(".*se=", "", l) else sub(".*se_t=", "", l)
    out[[test]] <- c(out[[test]], list(list(tag = tag, v = num(v))))
  }
  out
}
same1 <- function(a, b) {
  if (!identical(a$tag, b$tag) || length(a$v) != length(b$v)) return(FALSE)
  if (!identical(is.na(a$v), is.na(b$v))) return(FALSE)
  ok <- !is.na(a$v)
  !any(ok) || max(abs(a$v[ok] - b$v[ok]) / pmax(abs(a$v[ok]), 1e-300)) <= 1e-8
}
sameseq <- function(a, b) {
  length(a) == length(b) && all(vapply(seq_along(a), function(i) same1(a[[i]], b[[i]]), NA))
}
files <- sub("[.]trace:.*", "", sub("^frmtmb--", "",
  grep("keyed records differ", readLines(file.path(root, "rel083-log/se-cmp.txt")), value = TRUE)))
for (f in files) {
  tn <- file.path(root, "rel083-log/se-r22", paste0("frmtmb--", f, ".trace"))
  tb <- file.path(root, "rel083-log/se-r21", paste0("frmtmb--", f, ".trace"))
  N <- bytest(tn); B <- bytest(tb)
  L <- lapply(lanes, function(l) bytest(file.path(root, "reviews/rel083-review-out",
                                                  paste0("se-", l), paste0("frmtmb--", f, ".trace"))))
  names(L) <- lanes
  cat("FILE", f, "\n")
  for (t in union(names(N), names(B))) {
    if (!is.null(N[[t]]) && !is.null(B[[t]]) && sameseq(N[[t]], B[[t]])) next
    if (is.null(B[[t]])) { cat("  NEWTEST", t, "\n"); next }
    if (is.null(N[[t]])) { cat("  GONETEST", t, "\n"); next }
    asnew <- names(L)[vapply(L, function(x) !is.null(x[[t]]) && sameseq(x[[t]], N[[t]]), NA)]
    asold <- names(L)[vapply(L, function(x) !is.null(x[[t]]) && sameseq(x[[t]], B[[t]]), NA)]
    absent <- names(L)[vapply(L, function(x) is.null(x[[t]]), NA)]
    other <- setdiff(lanes, c(asnew, asold, absent))
    cat("  TEST", t, "| n", length(N[[t]]), "base", length(B[[t]]),
        "| as merged:", paste(asnew, collapse = " "), "| as 0.82.0:",
        length(asold), "| other:", paste(other, collapse = " "),
        if (length(absent)) paste("| absent:", paste(absent, collapse = " ")), "\n")
  }
}
