p <- "R/autocor.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("#' @rdname autocor-terms
#' @export
unstr <- function(time, gr) {
  ac_term_object(\"unstr\", sys.call())
}",
"#' @rdname autocor-terms
#' @export
unstr <- function(time, gr) {
  ac_term_object(\"unstr\", sys.call())
}

#' @rdname autocor-terms
#' @export
fcor <- function(M) {
  ac_term_object(\"fcor\", sys.call())
}")
rep1("#' @param autocor A one-sided formula of autocorrelation terms.",
"#' @param M For `fcor()`, the name of the known residual correlation or
#'   covariance matrix over the observations, which belongs in `data2`.
#' @param autocor A one-sided formula of autocorrelation terms.")
rep1("#' @return `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` return an
#'   object of class `frmtmb_ac_term` holding the term as a call.",
"#' @return `ar()`, `ma()`, `arma()`, `cosy()`, `unstr()` and `fcor()`
#'   return an object of class `frmtmb_ac_term` holding the term as a
#'   call.")
rep1("#' `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` are read inside a
#' model formula,",
"#' `ar()`, `ma()`, `arma()`, `cosy()`, `unstr()` and `fcor()` are read
#' inside a model formula,")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
