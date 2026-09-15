root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

# nlme owns fixef, ranef and VarCorr and is Recommended; generics owns
# refit and depends on base packages alone.  Both move out of Suggests.
sub1("DESCRIPTION",
paste0("Imports:\n    graphics,\n    grDevices,\n    Matrix,\n",
       "    methods,\n    mgcv,\n"),
paste0("Imports:\n    generics,\n    graphics,\n    grDevices,\n",
       "    Matrix,\n    methods,\n    mgcv,\n    nlme,\n"))
sub1("DESCRIPTION", "    mvtnorm,\n    nlme,\n", "    mvtnorm,\n")
cat("DONE\n")
