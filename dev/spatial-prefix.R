# Prefix, in place, every call in an R file to a frmtmb function that is
# not exported with frmtmb:::, for a test helper copied from core's suite
# into frmtmb.spatial's, which sees frmtmb's exports only. car_* and
# spde_* are frmtmb.spatial's own and are left alone.
#   Rscript dev/spatial-prefix.R <file>
f <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
x <- readLines(f, warn = FALSE)
ns <- asNamespace("frmtmb")
exn <- getNamespaceExports("frmtmb")
pd <- utils::getParseData(parse(text = x, keep.source = TRUE))
defined <- pd$text[pd$token == "SYMBOL" &
                     c(pd$token[-1], "") %in% c("LEFT_ASSIGN", "EQ_ASSIGN")]
calls <- unique(pd$text[pd$token == "SYMBOL_FUNCTION_CALL"])
internal <- calls[vapply(calls, function(g) {
  exists(g, envir = ns, inherits = FALSE) && !g %in% exn &&
    !g %in% defined && !grepl("^(car|spde)_", g)
}, NA)]
for (g in internal) {
  x <- gsub(paste0("(?<![:A-Za-z0-9_.$])", gsub(".", "[.]", g, fixed = TRUE),
                   "\\("), paste0("frmtmb:::", g, "("), x, perl = TRUE)
}
writeLines(x, f)
cat(f, "prefixed:", paste(internal, collapse = ", "), "\n")
