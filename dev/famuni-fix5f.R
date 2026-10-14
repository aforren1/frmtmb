p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("#' brms's `scale_xi()`: the GEV shape from an unbounded `tmp_xi`,
#' `lb + (ub - lb) plogis(tmp_xi)`, where `lb` and `ub` are the two values
#' `-1 / min(x)` and `-1 / max(x)` for `x = (y - mu) / sigma`, in order:
#' the interval of `xi` that keeps every observed row inside the support.
#'
#' @noRd
gev_scale_xi <- function(tmp, y, mu, sigma) {
  x <- (y - mu) / sigma
  b1 <- -1 / ad_tree_extreme(x, ad_min)
  b2 <- -1 / ad_tree_extreme(x, ad_max)
  lb <- ad_min(b1, b2)
  ub <- ad_max(b1, b2)
  lb + (ub - lb) / (1 + exp(-tmp))
}",
"#' brms's `scale_xi()`: the GEV shape from an unbounded `tmp_xi`,
#' `lb + (ub - lb) plogis(tmp_xi)`, where `lb` and `ub` are the two values
#' `-1 / min(x)` and `-1 / max(x)` for `x = (y - mu) / sigma`, in order:
#' the interval of `xi` that keeps every observed row inside the support.
#'
#' With one change: `lb` is held at -1 or above. brms's interval reaches
#' below -1, where the likelihood is unbounded; frmtmb's `xi` is above -1
#' on every route, and the fit's `log1p` link, which the draws are carried
#' back to, cannot represent a value below it. On the review's data at
#' xi = -0.85 (seed 1402, n = 40) brms's 95 percent interval reaches
#' -1.238; a first version without the hold gave a NaN `xi` for every draw
#' below -1 (dev/famuni-out/tmpxi-sample.txt). `ub` is held at `lb` or
#' above for the case where every row lies above `mu`, where the two
#' values are both below 0.
#'
#' @noRd
gev_scale_xi <- function(tmp, y, mu, sigma) {
  x <- (y - mu) / sigma
  b1 <- -1 / ad_tree_extreme(x, ad_min)
  b2 <- -1 / ad_tree_extreme(x, ad_max)
  lb <- ad_max(ad_min(b1, b2), -1)
  ub <- ad_max(ad_max(b1, b2), lb)
  lb + (ub - lb) / (1 + exp(-tmp))
}")
writeLines(strsplit(x, "\n")[[1]], p)
p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
r("      add(natural_dpar_prior(\"normal(0, 2.5)\", \"xi\", resp = rs))",
  "      add(set_prior(\"normal(0, 2.5)\", class = \"xi\", resp = rs))")
writeLines(strsplit(x, "\n")[[1]], p)
