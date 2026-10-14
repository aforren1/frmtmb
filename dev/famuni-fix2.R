p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  stopifnot(grepl(old, x, fixed = TRUE))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("#' The log of the Frechet scale, `log(mu) - lgamma(1 - 1 / nu)`, and a
#' term that is `NaN` where `nu <= 1` under a link that allows it. There
#' `gamma(1 - 1 / nu)` is negative or infinite, and `lgamma()` would
#' return the log of its absolute value, a finite and wrong scale.
#'
#' @noRd
frechet_log_scale <- function(dpars, lk, lk_nu) {
  nu <- dpars[[\"nu\"]]
  out <- dpar_log(dpars, \"mu\", lk) - lgamma(1 - 1 / nu)
  if (identical(lk_nu$name, \"logm1\")) return(out)
  out + 0 * log(nu - 1)
}",
"#' 0 where `x > 0` and `NaN` where `x <= 0`, with derivatives 0 where it\n#' is 0, on the RTMB tape as on doubles: a term that takes a parameter\n#' out of the likelihood where the density does not exist.\n#'\n#' With `u = x + |x|` it is `u / u - 1`. It must not have a constant 0\n#' as an operand: RTMB folds `0 * x` and `0 / x` to the constant 0 when\n#' it records the tape, so `0 * log(x)` is NaN on doubles and 0 in every\n#' fit (dev/rtmb-pitfalls.md, item 22). It raises no warning on doubles,\n#' where `log(x) - log(x)` would.\n#'\n#' @noRd\nnan_unless_pos <- function(x) {\n  u <- x + abs(x)\n  u / u - 1\n}\n\n#' The log of the Frechet scale, `log(mu) - lgamma(1 - 1 / nu)`, and a
#' term that is `NaN` where `nu <= 1` under a link that allows it. There
#' `gamma(1 - 1 / nu)` is negative or infinite, and `lgamma()` would
#' return the log of its absolute value, a finite and wrong scale.
#'
#' @noRd
frechet_log_scale <- function(dpars, lk, lk_nu) {
  nu <- dpars[[\"nu\"]]
  out <- dpar_log(dpars, \"mu\", lk) - lgamma(1 - 1 / nu)
  if (identical(lk_nu$name, \"logm1\")) return(out)
  out + nan_unless_pos(nu - 1)
}")
# the GEV: xi <= -1 has no likelihood under the identity link
r("  lk_xi <- dpar_link(link_xi, \"xi\", \"gen_extreme_value\", dpar_links_xi)
",
"  lk_xi <- dpar_link(link_xi, \"xi\", \"gen_extreme_value\", dpar_links_xi)
  # log1p holds xi above -1 by construction. Under the identity link
  # nothing does, and below -1 the likelihood is unbounded, so the soft
  # support, which caps an edge row at t = 1e-15, would give a finite
  # and wrong value there: seeds 621 and 622 ran to xi = -1.007 with a
  # row outside the support (dev/reviews/2026-10-13-famuni.md, fix 3)
  xi_guard <- if (identical(lk_xi$name, \"log1p\")) {
    function(dpars) 0
  } else {
    function(dpars) nan_unless_pos(1 + dpars[[\"xi\"]])
  }
")
r("      ll <- -dpar_log(dpars, \"sigma\", lk_sigma) - g[[\"l\"]] - g[[\"h\"]] -
        exp_neg_capped(g[[\"h\"]]) - 1e8 * g[[\"out\"]]",
"      ll <- -dpar_log(dpars, \"sigma\", lk_sigma) - g[[\"l\"]] - g[[\"h\"]] -
        exp_neg_capped(g[[\"h\"]]) - 1e8 * g[[\"out\"]] + xi_guard(dpars)")
r("      exp(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]]))",
  "      exp(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]])) + xi_guard(dpars)")
r("      log(-expm1(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]])))",
  "      log(-expm1(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]]))) +\n        xi_guard(dpars)")
writeLines(strsplit(x, "\n")[[1]], p)
