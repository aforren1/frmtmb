# Each measured pair of an outdir against the survey's record
# (dev/equiv-data/results.tsv): the measured verdict and level, and the
# numbers they rest on. Prints a table; writes nothing.
#
#   Rscript dev/translator-compare.R <outdir> [<id> ...]
#
# With no ids, every <id>.tsv in the outdir is compared.
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
ids <- args[-1]
rd <- function(p) read.delim(p, stringsAsFactors = FALSE,
                             na.strings = c("NA"))
rec <- rd("dev/equiv-data/results.tsv")
if (!length(ids)) {
  ids <- sub("[.]tsv$", "", list.files(outdir, pattern = "^p[0-9]+[.]tsv$"))
}
g <- function(d, nm) if (is.null(d[[nm]]) || is.na(d[[nm]][1])) "" else
  as.character(d[[nm]][1])
num <- function(d, nm) {
  v <- suppressWarnings(as.numeric(g(d, nm)))
  if (is.na(v)) "" else format(signif(v, 3))
}
rows <- lapply(ids, function(id) {
  p <- file.path(outdir, paste0(id, ".tsv"))
  r <- rec[rec$id == id, , drop = FALSE]
  if (!file.exists(p)) {
    return(data.frame(id = id, rec_verdict = g(r, "verdict_measured"),
                      rec_level = g(r, "level"), new_verdict = "NO RESULT",
                      new_level = "", same = FALSE, l1_const = "",
                      l1_grad = "", l1_extra = "", l1_marg = "", l1_dec = "",
                      l4_z = "", note = ""))
  }
  n <- rd(p)
  # the record keeps the measured verdict beside an adjudicated one
  rv <- g(r, "verdict_measured")
  rl <- g(r, "level")
  if (nzchar(g(r, "adjudicated"))) {
    rl <- paste0(rl, " (adj)")
  }
  data.frame(
    id = id, pair = paste(g(n, "feature_a"), "x", g(n, "feature_b")),
    rec_verdict = rv, rec_level = rl,
    new_verdict = g(n, "verdict"), new_level = g(n, "level"),
    same = identical(rv, g(n, "verdict")) &&
      identical(sub(" [(]adj[)]$", "", rl), g(n, "level")),
    rec_l1_const = num(r, "l1_const"), l1_const = num(n, "l1_const"),
    l1_grad = num(n, "l1_grad"), l1_extra = g(n, "l1_extra"),
    l1_marg = g(n, "l1_marg"), l1_dec = g(n, "l1_dec"),
    l4_z = g(n, "l4_z"),
    note = substr(paste(g(n, "l1_error"), g(n, "l4_error")), 1, 90),
    stringsAsFactors = FALSE)
})
allc <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(r) {
  for (m in setdiff(allc, names(r))) r[[m]] <- ""
  r[, allc]
})
tab <- do.call(rbind, rows)
op <- options(width = 400)
print(tab[, setdiff(names(tab), "note")], row.names = FALSE, right = FALSE)
nt <- tab[nzchar(trimws(tab$note)), c("id", "note")]
if (nrow(nt)) {
  cat("\nnotes (L1 or L4 errors):\n")
  for (i in seq_len(nrow(nt))) cat(" ", nt$id[i], nt$note[i], "\n")
}
options(op)
cat("\nsame verdict and level as the record:", sum(tab$same), "of",
    nrow(tab), "\n")
# the numbers too: a pair can keep its verdict while its constant moves
drift <- function(col) {
  d <- vapply(ids, function(id) {
    p <- file.path(outdir, paste0(id, ".tsv"))
    r <- rec[rec$id == id, , drop = FALSE]
    if (!file.exists(p) || !nrow(r)) return(NA_real_)
    n <- rd(p)
    a <- suppressWarnings(as.numeric(n[[col]][1]))
    b <- suppressWarnings(as.numeric(r[[col]][1]))
    if (is.null(a) || is.null(b) || length(a) != 1L || length(b) != 1L) {
      return(NA_real_)
    }
    abs(a - b)
  }, 0)
  cat(sprintf("%s: in both %d, bitwise equal %d, max |new - record| %s\n",
              col, sum(!is.na(d)), sum(d == 0, na.rm = TRUE),
              if (any(!is.na(d))) format(max(d, na.rm = TRUE)) else "NA"))
}
for (col in c("l1_const", "l1_grad", "l2_metric")) drift(col)
