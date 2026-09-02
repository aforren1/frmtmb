# G5.2a bookkeeping: every condition message raised in R/ must be
# unique. This walks the parsed ASTs, joins the literal string pieces of
# each stop()/warning() call (interpolated values are dropped, since two
# messages that differ only by an interpolated value are still the same
# sentence), and reports duplicates.
#
# Run: Rscript dev/tre/count-messages.R [dir]
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "R"

lits <- function(e) {
  if (is.character(e)) return(e)
  if (!is.call(e)) return(character(0))
  unlist(lapply(as.list(e)[-1], lits))
}

msgs <- list()
walk <- function(e, file) {
  if (!is.call(e)) return(invisible())
  h <- e[[1]]
  if (is.name(h) && as.character(h) %in% c("stop", "warning")) {
    a <- as.list(e)[-1]
    nm <- names(a) %||% rep("", length(a))
    a <- a[!nm %in% c("call.", "domain", "immediate.")]
    txt <- paste(unlist(lapply(a, lits)), collapse = "")
    if (nzchar(trimws(txt))) {
      msgs[[length(msgs) + 1L]] <<- list(file = file, msg = txt)
    }
  }
  for (i in seq_along(e)) {
    if (is.call(e[[i]])) walk(e[[i]], file)
  }
  invisible()
}
`%||%` <- function(a, b) if (is.null(a)) b else a

for (f in list.files(dir, pattern = "\\.R$", full.names = TRUE)) {
  for (e in parse(f)) walk(e, basename(f))
}

tab <- data.frame(file = vapply(msgs, `[[`, "", "file"),
                  msg = vapply(msgs, `[[`, "", "msg"),
                  stringsAsFactors = FALSE)
cat("condition-raising calls: ", nrow(tab), "\n", sep = "")
cat("distinct messages:       ", length(unique(tab$msg)), "\n", sep = "")
dup <- tab$msg[duplicated(tab$msg)]
if (length(dup)) {
  cat("\nDUPLICATES (", length(unique(dup)), "):\n", sep = "")
  for (m in unique(dup)) {
    cat("  [", paste(tab$file[tab$msg == m], collapse = ", "), "] ",
        substr(m, 1, 90), "\n", sep = "")
  }
} else {
  cat("no duplicates\n")
}
