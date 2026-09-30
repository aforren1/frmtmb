# Lane sampfix nits round: a probe that fails at both fills re-raises.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/R/methods-draws.R"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, ":\n", old)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("#' `0` read them. `NA` if both fills fail, which is an error about
#' something else.",
"#' `0` read them. `NA` if both fills fail, which is an error about
#' something else; the error rides along as the `error` attribute.")
rep1("  if (inherits(a, \"error\") && inherits(b, \"error\")) return(NA)",
"  if (inherits(a, \"error\") && inherits(b, \"error\")) {
    return(structure(NA, error = a))
  }")
rep1("draws_laplace_probe <- function(x, what, at, at_na = NULL) {
  if (!draws_is_laplace(x)) return(invisible(NULL))
  if (!isTRUE(draws_probe_reads(at))) return(invisible(NULL))
  draws_laplace_refuse(x, what, at_na)
}",
"draws_laplace_probe <- function(x, what, at, at_na = NULL) {
  if (!draws_is_laplace(x)) return(invisible(NULL))
  r <- draws_probe_reads(at)
  # an error at both fills is the call's own, and the first draw of the
  # computation would raise it anyway; raising it here keeps a probe that
  # fails on its own account from waving the call through
  if (is.na(r)) stop(attr(r, \"error\"))
  if (!r) return(invisible(NULL))
  draws_laplace_refuse(x, what, at_na)
}")
writeLines(s, f)
cat("ok\n")
