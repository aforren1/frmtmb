# Derive, in its own process off the BASE install, the set of S3
# generics frmtmb.sample DEFINES. Independent of the lane's table.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- loadNamespace("frmtmb.sample")
info <- parseNamespaceFile("frmtmb.sample",
                           "C:/Users/adf44/source/r/rellib-r3")
exp_names <- sort(unique(info$exports))
out <- character()
bodies <- list()
frms <- list()
for (nm in exp_names) {
  f <- tryCatch(get(nm, envir = ns, inherits = FALSE),
                error = function(e) NULL)
  if (!is.function(f)) next
  if (!identical(environment(f), ns)) next
  b <- paste(deparse(body(f)), collapse = " ")
  if (grepl("UseMethod", b, fixed = TRUE)) {
    out <- c(out, nm)
    bodies[[nm]] <- b
    frms[[nm]] <- paste(deparse(args(f)), collapse = " ")
  }
}
# also: every export that is a function, and its origin, so the
# re-export set can be separated from the defined set
origin <- vapply(exp_names, function(nm) {
  f <- tryCatch(get(nm, envir = ns, inherits = FALSE),
                error = function(e) NULL)
  if (!is.function(f)) return(NA_character_)
  e <- environment(f)
  if (is.null(e)) "base" else environmentName(e)
}, "")
dir.create("C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/sgrev-out",
           showWarnings = FALSE, recursive = TRUE)
saveRDS(list(names = out, bodies = bodies, formals = frms,
             exports = exp_names, origin = origin),
        "C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/sgrev-out/names.rds")
cat("defined generics:", length(out), "\n")
cat(paste(out, collapse = " "), "\n")
cat("\nbodies that are NOT a bare UseMethod:\n")
bare <- vapply(out, function(nm) {
  b <- gsub("[[:space:]]+", "", bodies[[nm]])
  grepl("^\\{?UseMethod\\(\"", b)
}, NA)
cat(paste(out[!bare], collapse = " "), "\n")
for (nm in out[!bare]) cat("  ", nm, ": ", bodies[[nm]], "\n", sep = "")
