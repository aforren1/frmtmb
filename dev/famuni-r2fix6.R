p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("#' it a supremum. So a scalar `xi` is judged by the sign of
#' `d logLik / d xi` at the stop, from the objective's gradient carried
#' to the `xi` scale: negative, the likelihood is still rising toward -1
#' and the warning says so; otherwise the fit is at a maximum inside the
#' band and gets a message instead. A predicted `xi` has no one
#' direction to read, and gets a warning that names the distance only.",
"#' it a supremum. So a scalar `xi` is judged by `d logLik / d xi` at the
#' stop, from the objective's gradient carried to the `xi` scale. Where
#' it is not near 0 the edge, not a maximum, stopped the fit, and the
#' warning says so: negative, the likelihood still rises toward -1;
#' large and positive, the fit sits in the corner where moving `xi` up
#' would push the largest row out of the support, which the identity
#' link reaches (3 of the 720 fits, at xi + 1 below 1e-15, d logLik /
#' d xi 154 to 238). Near 0 the fit is at a maximum inside the band and
#' gets a message instead. On the 720 fits (dev/famuni-r2-xiedge.R) the
#' edge fits have `|d logLik / d xi|` of 0.88 or more and the one
#' interior maximum 3.3e-4; 0.1 separates them with room on both sides.
#' A predicted `xi` has no one direction to read, and gets a warning
#' that names the distance only.")
r("  } else if (dll < 0) {
    frm_warning(\"gen_extreme_value: xi \", of, \"ran to \", at, \", the edge \",
                \"at -1, with the likelihood still rising toward it (d \",
                \"logLik / d xi = \", format(dll, digits = 3), \"). Below -1 \",",
"  } else if (abs(dll) > 0.1 || dll < 0) {
    frm_warning(\"gen_extreme_value: xi \", of, \"ran to \", at, \", the edge \",
                \"at -1, and the likelihood is not stationary in xi there \",
                \"(d logLik / d xi = \", format(dll, digits = 3), \"): the \",
                \"edge, not a maximum, stopped the fit. Below -1 \",")
writeLines(strsplit(x, "\n")[[1]], p)
