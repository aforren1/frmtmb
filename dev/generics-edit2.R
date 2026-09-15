root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  txt <- sub(old, new, txt, fixed = TRUE)
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

sub1("R/methods-fit.R",
paste0("#' Extract random-effect covariance matrices\n",
       "#' @param x A `frmtmb_fit`.\n",
       "#' @param ... Unused.\n"),
paste0("#' Extract random-effect covariance matrices\n",
       "#' @param x A `frmtmb_fit`.\n",
       "#' @param sigma Ignored. It is carried by nlme's generic, which\n",
       "#'   frmtmb now shares rather than shadows, for the models that\n",
       "#'   scale a covariance by a residual standard deviation.\n",
       "#' @param ... Unused.\n"))
cat("DONE\n")
