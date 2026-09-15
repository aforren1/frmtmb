root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L)
    stop("many in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

nl <- "\\n"
old <- paste0(
'cat(sprintf("%-20s %-12s %-12s %-10s', nl, '", "generic", "generic-from",\n',
'            "method-from", "matched"))\n',
'for (g in gens) {\n',
'  f <- tryCatch(get(g), error = function(e) NULL)\n',
'  r <- resolve(g, cls)\n',
'  cat(sprintf("%-20s %-12s %-12s %-10s', nl, '", g, where_from(f),\n',
'              r[["pkg"]], r[["how"]]))\n',
'}')
new <- paste0(
'cat(sprintf("%-20s %-12s %-12s %-10s %-6s', nl, '", "generic",\n',
'            "generic-from", "method-from", "matched", "getS3"))\n',
'for (g in gens) {\n',
'  f <- tryCatch(get(g), error = function(e) NULL)\n',
'  r <- resolve(g, cls)\n',
'  gs <- !is.null(tryCatch(getS3method(g, cls, optional = TRUE),\n',
'                          error = function(e) NULL))\n',
'  cat(sprintf("%-20s %-12s %-12s %-10s %-6s', nl, '", g, where_from(f),\n',
'              r[["pkg"]], r[["how"]], gs))\n',
'}')
sub1("dev/generics-check.R", old, new)
cat("DONE\n")
