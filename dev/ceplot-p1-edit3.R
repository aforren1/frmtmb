# Lane ceplot punch 1 (m7): frm_check_dots() takes `.hidden`, names a
# call accepts without listing them among its arguments (a deprecated
# alias such as plot()'s do_plot).
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/utils.R"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(s, old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("expected one match, got ", n, ": ", substr(old, 1, 60))
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, "#'   retired spelling is still refused by a method that forwards.
#' @noRd
frm_check_dots <- function(..., .unsupported = NULL, .allow = NULL) {",
"#'   retired spelling is still refused by a method that forwards.
#' @param .hidden Character vector of names accepted like `.allow` but
#'   not listed back in the refusal's \"It takes:\", for a deprecated
#'   alias that the caller should not be pointed to.
#' @noRd
frm_check_dots <- function(..., .unsupported = NULL, .allow = NULL,
                           .hidden = NULL) {")
s <- rep1(s, "  bad <- setdiff(bad, .allow)
  if (!length(bad)) return(invisible(NULL))",
"  bad <- setdiff(bad, c(.allow, .hidden))
  if (!length(bad)) return(invisible(NULL))")
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
