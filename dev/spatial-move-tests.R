# Move test_that() blocks (and named helper definitions) from a core
# test file to a frmtmb.spatial test file, verbatim, with the comment
# lines directly above each. Calls to frmtmb functions that are not
# exported get a frmtmb::: prefix, because the moved file runs in
# frmtmb.spatial's test environment, which sees frmtmb's exports only.
#
#   Rscript dev/spatial-move-tests.R <from> <to> <pattern> [helpers]
#
# <pattern> is a regular expression on the test_that() description;
# [helpers] a comma-separated list of top-level object names to move
# (they are COPIED, not removed, when a block left behind still uses
# them; say "copy:" before the list to copy). The target is appended
# to. Prints what moved.
a <- commandArgs(trailingOnly = TRUE)
from <- a[1]; to <- a[2]; pat <- a[3]
hl <- if (length(a) > 3) a[4] else ""
copy_helpers <- startsWith(hl, "copy:")
helpers <- strsplit(sub("^copy:", "", hl), ",", fixed = TRUE)[[1]]
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
x <- readLines(from, warn = FALSE)
ex <- parse(text = x, keep.source = TRUE)
sr <- attr(ex, "srcref")
take <- list(); keep_lines <- rep(TRUE, length(x))
prev_end <- 0L
for (i in seq_along(ex)) {
  e <- ex[[i]]
  s <- sr[[i]][1]; f <- sr[[i]][3]
  k <- s - 1L
  while (k > prev_end && grepl("^\\s*#", x[k])) k <- k - 1L
  prev_end <- f
  is_test <- is.call(e) && identical(e[[1]], as.name("test_that"))
  nm <- if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
            is.name(e[[2]])) as.character(e[[2]]) else NA
  hit_test <- is_test && grepl(pat, as.character(e[[2]]))
  hit_help <- !is.na(nm) && nm %in% helpers
  if (hit_test || hit_help) {
    take[[length(take) + 1L]] <- x[(k + 1L):f]
    cat(if (hit_test) "TEST" else "HELPER", ":",
        if (hit_test) as.character(e[[2]]) else nm, "\n")
    if (hit_test || !copy_helpers) keep_lines[(k + 1L):f] <- FALSE
  }
}
moved <- unlist(lapply(take, function(t) c(t, "")))
# prefix frmtmb internals called in the moved code
ns <- asNamespace("frmtmb")
ex_names <- getNamespaceExports("frmtmb")
pd <- utils::getParseData(parse(text = moved, keep.source = TRUE))
calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", ]
internal <- unique(calls$text[vapply(calls$text, function(f) {
  exists(f, envir = ns, inherits = FALSE) && !f %in% ex_names &&
    !grepl("^(car|spde)_", f)
}, NA)])
# a call already written frmtmb:::f keeps its prefix
for (f in internal) {
  moved <- gsub(paste0("(?<![:A-Za-z0-9_.$])", gsub(".", "[.]", f,
                                                       fixed = TRUE),
                       "\\("),
                paste0("frmtmb:::", f, "("), moved, perl = TRUE)
}
cat("prefixed:", paste(internal, collapse = ", "), "\n")
cat(moved, file = to, sep = "\n", append = TRUE)
writeLines(x[keep_lines], from)
