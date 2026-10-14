p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("#' `link_nu = \"identity\"` a `nu` at or below 1 has no mean, and the
#' density is `NaN` there, as brms's is.",
"#' `link_nu = \"identity\"` a `nu` at or below 1 has no mean, and the
#' density is `NaN` there, as brms's is, in the fitted objective as well
#' as on doubles.")
r("#' that turns it back, not a wall. The start, `xi = 0`, is inside the
#' support for every data set, and a fit with one component ends with
#' every row inside it, the region brms's `scale_xi()` maps its `tmp_xi`
#' onto.",
"#' that turns it back, not a wall. The start, `xi = 0`, is inside the
#' support for every data set. For `xi > -1` the density goes to 0 at
#' an end of the support, so a maximum has every row inside it. That is
#' the region brms's `scale_xi()` maps its `tmp_xi` onto.")
r("#' The default link of `xi` is `log1p`, `log(1 + xi)`, brms's link
#' for a predicted `xi`, which keeps `xi > -1`; `link_xi = \"identity\"`
#' is the other choice brms allows. brms parameterizes an `xi` without
#' a formula differently, as an unbounded `tmp_xi` scaled into the
#' support the current `mu` and `sigma` leave, and reports `xi` from it
#' after the fit. That is a change of variable, so the maximum of the
#' likelihood is the same point. It is not used here because its range
#' reaches below `xi = -1`, where the likelihood is unbounded and has no
#' maximum, so frmtmb holds an `xi` without a formula to the `log1p`
#' link as well. For `-1 < xi < -0.5` the maximum exists but is not
#' regular (Smith, 1985), so its standard errors do not have their usual
#' meaning.",
"#' The default link of `xi` is `log1p`, `log(1 + xi)`, brms's link
#' for a predicted `xi`, which keeps `xi > -1`. `link_xi = \"identity\"`
#' is the other choice brms allows. Under it, `xi <= -1` has no
#' likelihood, so every link holds `xi` above -1. Below -1 the
#' likelihood is unbounded and has no maximum.
#'
#' brms parameterizes an `xi` without a formula differently. It samples
#' an unbounded `tmp_xi` and scales it into the support the current `mu`
#' and `sigma` leave, and it reports `xi` from it after the fit. That is
#' a change of variable, so the maximum of the likelihood is the same
#' point. A fit does not use it, because its range reaches below
#' `xi = -1`. [frmtmb.sample::frm_sample()] does, because there the
#' support's edge is a wall the sampler meets; it holds the range above
#' -1 (see [frmtmb-gev-xi-chart]).
#'
#' A fit whose `xi` stops within 1e-3 of -1 warns: the likelihood rises
#' toward the edge, so its log-likelihood is a supremum, not a maximum.
#' For `-1 < xi < -0.5` the maximum exists but is not regular (Smith,
#' 1985), so its standard errors do not have their usual meaning; no
#' message says so.")
r("#' `cens()` follows frmtmb's discrete censoring convention (see
#' [frmtmb_family()], Censoring a discrete response), and `trunc(lb = )`
#' reads its bound as inclusive, as brms does.",
"#' `cens()` follows frmtmb's discrete censoring convention (see
#' [frmtmb_family()], Censoring a discrete response): a right-censored
#' count at `k` is `P(Y >= k)`. brms reads it as `P(Y > k)`, so right
#' and interval censoring differ from brms; left censoring agrees.
#' `trunc(lb = )` reads its bound as inclusive, as brms does.")
r("    # an ordinary gradient that points back in. Inside the support the
    # term is exactly 0.",
"    # an ordinary gradient that points back in. It is exactly 0 where
    # t > 1e-15, which the clamp in gev_parts() starts at.")
writeLines(strsplit(x, "\n")[[1]], p)
