# Lane sampfix nits round: restore is_na_re_form(), dropped in the rewrite
# of the probe block and caught by R CMD check.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/R/methods-draws.R"
s <- gsub("\r", "", readLines(f))
i <- which(s == "#' A check on every draw's result, after `draws_laplace_probe()` let the")
stopifnot(length(i) == 1L)
add <- c(
"#' Whether `re_formula` is the `NA` that drops every group-level term.",
"#'",
"#' @noRd",
"is_na_re_form <- function(re_form) {",
"  !is.null(re_form) && length(re_form) == 1L && !is.language(re_form) &&",
"    is.na(re_form)",
"}",
"")
writeLines(c(s[seq_len(i - 1L)], add, s[i:length(s)]), f)
cat("ok\n")
