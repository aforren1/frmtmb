# Short check, item 7: static checks on the trial merge. Top-level names defined twice in
# R/ (the later file in collation wins without a word), and whether frmtmb.sample's calls
# into core all resolve to exports of the merged core.
M <- "C:/Users/adf44/source/r/rev-predread-src/merge"
fs <- list.files(file.path(M, "R"), pattern = "[.]R$", full.names = TRUE)
defs <- list()
for (f in fs) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && (identical(e[[1L]], as.name("<-")) || identical(e[[1L]], as.name("="))) && is.name(e[[2L]])) {
      nm <- as.character(e[[2L]])
      defs[[nm]] <- c(defs[[nm]], basename(f))
    }
  }
}
dup <- defs[lengths(defs) > 1L]
cat("top-level names defined more than once:", length(dup), "\n")
for (n in names(dup)) cat(" ", n, ":", paste(dup[[n]], collapse = ", "), "\n")
b <- function(path, s) {
  x <- list.files(path, pattern = "[.]R$", full.names = TRUE)
  unique(unlist(lapply(x, function(f) {
    ex <- parse(f, keep.source = FALSE)
    unlist(lapply(ex, function(e) if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]])) as.character(e[[2L]])))
  })))
}
ns <- readLines(file.path(M, "NAMESPACE"))
txt <- paste(ns, collapse = " ")
cat("markers in NAMESPACE:", any(grepl("^<<<<<<<|^>>>>>>>", ns)), "\n")
.libPaths(c("C:/Users/adf44/source/r/rev-predread-lib/merge", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
ex <- getNamespaceExports("frmtmb")
sample_defs <- b(file.path(M, "extensions/frmtmb.sample/R"))
calls <- unique(unlist(lapply(list.files(file.path(M, "extensions/frmtmb.sample/R"), pattern = "[.]R$", full.names = TRUE), function(f) {
  all.names(parse(f, keep.source = FALSE), functions = TRUE)
})))
core <- names(defs)
need <- setdiff(intersect(calls, core), c(sample_defs, ex))
cat("names frmtmb.sample uses that core defines and does not export:", length(need), "\n")
print(utils::head(need, 40))
for (ln in c("predread", "latres", "nlfix")) {
  nsf <- readLines(sprintf("C:/Users/adf44/source/r/frmtmb-wt-%s/NAMESPACE", ln))
  cat(ln, "NAMESPACE lines", length(nsf), "\n")
}
cat("merged export check, six of predread:", all(c("ce_na_vars", "ce_me_grouped", "autocor_for_newdata", "update_merge_prior", "ce_by_response", "resp_refuse_single") %in% ex), "\n")
