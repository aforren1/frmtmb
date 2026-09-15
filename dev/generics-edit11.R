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

sub1("extensions/frmtmb.sample/R/methods-draws.R",
paste0("#' @param object A `frmtmb_draws`, or a matrix of draws\n",
       "#'   (variables in columns).\n",
       "#' @param probs Quantiles for `posterior_summary()`.\n"),
paste0("#' @param object A `frmtmb_draws`, or a matrix of draws\n",
       "#'   (variables in columns).\n",
       "#' @param x The same, for `posterior_summary()`, whose generic\n",
       "#'   is brms's and names its first argument `x`.\n",
       "#' @param probs Quantiles for `posterior_summary()`.\n"))
cat("DONE\n")
