# The R CMD check sections a NAMESPACE, a signature or an Rd change can
# break, run directly on the installed packages: R CMD check itself is
# held for the release. Each prints nothing when clean.
#   Rscript dev/samplegen-toolscheck.R <LIB>
av <- commandArgs(trailingOnly = TRUE)
.libPaths(c(av[1], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("frmtmb", "frmtmb.sample")) {
  cat("====", p, "from", find.package(p), "\n")
  for (fn in c("undoc", "codoc", "checkDocFiles", "checkS3methods",
               "checkReplaceFuns")) {
    r <- tryCatch(utils::capture.output(print(
      get(fn, asNamespace("tools"))(package = p))),
      error = function(e) paste("ERROR:", conditionMessage(e)))
    r <- r[nzchar(r)]
    cat(sprintf("-- %-16s %s\n", fn,
                if (length(r)) paste(r, collapse = "\n   ") else "clean"))
  }
}
