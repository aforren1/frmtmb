# Reviewer, lane wt-conditions (round 2): every frm_stop(), frm_warning()
# and frm_message() call site, with whether it keeps the call. default,
# for choosing the fidelity sample.
#   Rscript dev/conditions-rev-sitetable.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
dirs <- c(file.path(root, "R"), Sys.glob(file.path(root, "extensions/*/R")),
          file.path(root, "inst"))
rows <- list()
for (d in dirs) {
  for (f in list.files(d, "[.]R$", full.names = TRUE, recursive = TRUE)) {
    ex <- parse(f, keep.source = TRUE)
    pd <- getParseData(ex)
    hit <- pd$token == "SYMBOL_FUNCTION_CALL" &
      pd$text %in% c("frm_stop", "frm_warning", "frm_message")
    for (id in pd$parent[hit]) {
      callid <- pd$parent[pd$id == id]
      txt <- getParseText(pd, callid)
      e <- str2lang(txt)
      hascall <- "call." %in% names(e)
      rows[[length(rows) + 1L]] <- data.frame(
        pkg = if (basename(d) == "R") basename(dirname(d)) else "inst",
        file = sub(paste0(d, "/"), "", f, fixed = TRUE),
        line = pd$line1[pd$id == callid], fn = as.character(e[[1L]]),
        call_default = !hascall,
        call_false = hascall && identical(e[["call."]], FALSE),
        nargs = length(e) - 1L - hascall,
        fmt = grepl("sprintf|gettextf|paste|format[(]|deparse", txt))
    }
  }
}
x <- do.call(rbind, rows)
print(table(x$fn, ifelse(x$call_false, "call.=FALSE",
                         ifelse(x$call_default, "default", "other"))))
print(table(x$pkg, x$call_default))
saveRDS(x, file.path(root, "dev/conditions-rev-log/sitetable.rds"))
