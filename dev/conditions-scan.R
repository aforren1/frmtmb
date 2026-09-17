# Lane wt-conditions: census of condition-raising call sites, from the
# parse tree, before any rewrite. Run from the worktree root.
roots <- c("R", file.path(list.dirs("extensions", recursive = FALSE), "R"))
names_of_interest <- c("stop", "warning", "message", "stopifnot",
                       "simpleError", "simpleWarning", "simpleMessage",
                       "simpleCondition", "signalCondition",
                       "errorCondition", "warningCondition", "abort",
                       "warn", "inform", ".Defunct", ".Deprecated",
                       "packageStartupMessage", "tryCatch",
                       "withCallingHandlers", "invokeRestart")
rows <- list()
for (r in roots) {
  for (f in list.files(r, pattern = "[.][Rr]$", full.names = TRUE)) {
    pd <- getParseData(parse(f, keep.source = TRUE))
    fc <- pd[pd$token %in% c("SYMBOL_FUNCTION_CALL", "SYMBOL"), ]
    fc <- fc[fc$text %in% names_of_interest, ]
    if (nrow(fc)) {
      rows[[length(rows) + 1L]] <- data.frame(root = r, file = f,
        line = fc$line1, token = fc$token, text = fc$text)
    }
  }
}
d <- do.call(rbind, rows)
print(with(d, table(root, text, token)))
nonCall <- d[d$token == "SYMBOL", ]
if (nrow(nonCall)) print(nonCall)
