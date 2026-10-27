# Review: roxygenise copies of core, frmtmb.sample and frmtmb.spatial
# (dev/review-rel076-out/rox/) against rellib-r14 and compare the
# regenerated man/ and NAMESPACE with the release tree's.
#   Rscript dev/review-rel076-rox.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
rel <- "C:/Users/adf44/source/r/frmtmb-wt-rel076"
rox <- file.path(rel, "dev/review-rel076-out/rox")
src <- c(frmtmb = rel,
         frmtmb.sample = file.path(rel, "extensions/frmtmb.sample"),
         frmtmb.spatial = file.path(rel, "extensions/frmtmb.spatial"))
for (p in names(src)) {
  d <- file.path(rox, p)
  # the copy keeps the tree's Rd: roxygen resolves a link to a local
  # topic through the man/ it finds (a first run from an empty man/
  # linked 64 core topics to brms), and deletes an Rd of its own header
  # that it no longer makes, so a stale page still shows
  suppressMessages(roxygen2::roxygenise(d))
  a <- sort(list.files(file.path(src[[p]], "man"), "[.]Rd$"))
  b <- sort(list.files(file.path(d, "man"), "[.]Rd$"))
  cat("==", p, "tree Rd", length(a), "regenerated", length(b), "\n")
  cat("  only tree:", setdiff(a, b), "\n  only regenerated:", setdiff(b, a),
      "\n")
  rd <- function(f) gsub("\r", "", readLines(f, warn = FALSE))
  ndiff <- 0
  for (f in intersect(a, b)) {
    if (!identical(rd(file.path(src[[p]], "man", f)),
                   rd(file.path(d, "man", f)))) {
      ndiff <- ndiff + 1
      cat("  DIFFERS:", f, "\n")
    }
  }
  same_ns <- identical(rd(file.path(src[[p]], "NAMESPACE")),
                       rd(file.path(d, "NAMESPACE")))
  cat("  Rd differ:", ndiff, " NAMESPACE equal:", same_ns, "\n")
}
