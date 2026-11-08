# Show that the guard in tests/testthat/test-no-superassign.R fails on
# the code it guards against: run its two scanners on the base build
# (rellib-r18, frmtmb 0.79.1 at 3cd2a0c8) and on the base sources.
#
# Usage: Rscript dev/superassign-guard-base.R <lib> <srcdir-with-R>

a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1L], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- a[2L]

# the scanners exactly as the test defines them: evaluate the test
# file's top-level function definitions and nothing else
ex <- parse("tests/testthat/test-no-superassign.R", keep.source = FALSE)
for (e in ex) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
        is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
    eval(e, globalenv())
  }
}

ns <- asNamespace("frmtmb")
cat("lib:", dirname(getNamespaceInfo(ns, "path")), "version:",
    format(utils::packageVersion("frmtmb")), "\n")
hits_ns <- superassign_ns_scan(ns)
cat(sprintf("NS functions with <<-: %d\n", length(hits_ns)))
cat(paste0("  ", hits_ns, collapse = "\n"), "\n")

hits_src <- superassign_src_scan(file.path(src, "R"))
cat(sprintf("SRC <<- tokens in %s/R: %d in %d files\n", src,
            length(hits_src),
            length(unique(sub(":.*", "", trimws(hits_src))))))
for (d in list.dirs(file.path(src, "extensions"), recursive = FALSE)) {
  h <- superassign_src_scan(file.path(d, "R"))
  cat(sprintf("SRC %s: %d\n", basename(d), length(h)))
}
