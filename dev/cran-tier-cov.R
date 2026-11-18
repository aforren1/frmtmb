# Lane cran: what the proposed CRAN tier exercises, export by export.
# Joins the trace counts of dev/cran-tier-run.R (calls that REACHED each
# export while the tier ran) with the names each kept block and each
# example uses, and splits the exports into user-facing ones and the
# extension API (aliases of frmtmb-sampling-api, frmtmb-extension-api,
# frmtmb-shared-generics and \keyword{internal} topics).
#
#   Rscript dev/cran-tier-cov.R <worktree> > dev/cran-log/tier-cov.txt
local({
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- commandArgs(trailingOnly = TRUE)[1]
logd <- file.path(root, "dev", "cran-log")
ex <- sort(unique(parseNamespaceFile("frmtmb", dirname(root))$exports))
calls <- utils::read.delim(file.path(logd, "tierrun-trace-calls.tsv"))
reached <- setNames(calls$calls, calls$export)

rd <- tools::Rd_db(dir = root)
api_topics <- c("frmtmb-sampling-api.Rd", "frmtmb-extension-api.Rd",
                "frmtmb-shared-generics.Rd")
topic_of <- list()
api <- character()
ex_names <- character()
for (nm in names(rd)) {
  r <- rd[[nm]]
  al <- unlist(tools:::.Rd_get_metadata(r, "alias"))
  kw <- unlist(tools:::.Rd_get_metadata(r, "keyword"))
  if (nm %in% api_topics || "internal" %in% kw) api <- c(api, al)
  for (a in al) topic_of[[a]] <- c(topic_of[[a]], sub("[.]Rd$", "", nm))
  f <- tempfile(fileext = ".R")
  tools::Rd2ex(r, f, commentDontrun = TRUE, commentDonttest = FALSE)
  if (file.exists(f)) {
    ex_names <- c(ex_names, tryCatch(all.names(parse(f, keep.source = FALSE)),
                                     error = function(e) character()))
    unlink(f)
  }
}
ex_names <- unique(ex_names)
kind <- ifelse(ex %in% api, "api", "user")
tr <- ifelse(ex %in% names(reached), reached[ex], NA)
d <- data.frame(export = ex, kind = kind, reached = tr,
                example = ex %in% ex_names)
d$status <- ifelse(!is.na(d$reached) & d$reached > 0, "reached",
                   ifelse(is.na(d$reached), "untraced", "not reached"))
cat("== exports by kind and tier status (trace of the tier run)\n")
print(table(d$kind, d$status))
cat("\n== user-facing exports the tier does not reach\n")
u <- d[d$kind == "user" & d$status != "reached", ]
for (j in seq_len(nrow(u))) {
  cat("USER\t", u$export[j], "\t", u$status[j], "\texample=", u$example[j],
      "\ttopic=", paste(topic_of[[u$export[j]]], collapse = ","), "\n",
      sep = "")
}
cat("\n== extension-API exports the tier does not reach\n")
v <- d[d$kind == "api" & d$status != "reached", ]
cat(paste0(v$export, " (", v$status, ")"), sep = ", ", fill = 78)
utils::write.table(d, file.path(logd, "tier-cov.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
})
