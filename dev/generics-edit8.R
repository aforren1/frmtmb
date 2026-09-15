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

# predict() is the method whose default scale differs from brms's, so
# its page carries the pointer to the contract
sub1("R/predict.R",
paste0("#' head(predict(fit2, dpar = \"sigma\", type = \"response\"))\n",
       "#' @export\npredict.frmtmb_fit <- function(object,"),
paste0("#' head(predict(fit2, dpar = \"sigma\", type = \"response\"))\n",
       "#' @seealso [frmtmb-scales], which states which scale every\n",
       "#'   method reports. The default here is the LINK scale, where\n",
       "#'   brms's `predict()` gives the response scale.\n",
       "#' @export\npredict.frmtmb_fit <- function(object,"))

sub1("R/predict.R",
paste0("#' r <- residuals(fit, type = \"osa\")\n",
       "#' qqnorm(r); qqline(r)\n#' @export\n"),
paste0("#' r <- residuals(fit, type = \"osa\")\n",
       "#' qqnorm(r); qqline(r)\n",
       "#' @seealso [frmtmb-scales] for which scale each type is on.\n",
       "#' @export\n"))
cat("DONE\n")
