# Review: roxygenise a scratch copy of the tree (argument 1) and report
# whether NAMESPACE and man/ change, for core and the five changed
# extensions. The copy is outside the tree; the tree is not touched.
args <- commandArgs(TRUE)
cp <- args[[1]]
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
dirs <- c(cp, file.path(cp, "extensions", c("frmtmb.survival",
  "frmtmb.spatial", "frmtmb.sample", "frmtmb.latent", "frmtmb.eam",
  "frmtmb.coupling", "frmtmb.ode", "frmtmb.learn")))
snap <- function(d) {
  f <- c(file.path(d, "NAMESPACE"), list.files(file.path(d, "man"),
                                                full.names = TRUE))
  setNames(vapply(f, function(x) paste(sub("\r$", "", readLines(x, warn = FALSE)),
                                       collapse = "\n"), ""), basename(f))
}
for (d in dirs) {
  a <- snap(d)
  res <- tryCatch({
    suppressMessages(roxygen2::roxygenise(d))
    "ok"
  }, error = function(e) conditionMessage(e))
  b <- snap(d)
  ch <- union(setdiff(names(a), names(b)), setdiff(names(b), names(a)))
  both <- intersect(names(a), names(b))
  ch <- c(ch, both[a[both] != b[both]])
  cat(basename(d), ": roxygenise", res, "; files changed:", length(ch),
      if (length(ch)) paste(ch, collapse = " ") else "", "\n")
}
