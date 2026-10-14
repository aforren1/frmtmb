p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  stopifnot(grepl(old, x, fixed = TRUE))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("fam_gen_extreme_value <- function(link = \"identity\", link_sigma = \"log\",
                                  link_xi = \"log1p\") {
  lk <- mu_link(link, \"gen_extreme_value\")",
"fam_gen_extreme_value <- function(link = \"identity\", link_sigma = \"log\",
                                  link_xi = \"log1p\") {
  fam_gev_impl(link, link_sigma, link_xi, chart = \"xi\")
}

#' The GEV family in one of two charts for `xi`. `\"xi\"` is the family
#' `gen_extreme_value()` returns. `\"tmp\"` is brms's: the `xi` dpar holds
#' `tmp_xi` on an identity link, and the density maps it into the support
#' the current `mu` and `sigma` leave (`gev_scale_xi()`), REPORTing the
#' `xi` it used as `frmtmb_gev_xi`. Only the sampler tapes it
#' (`gev_xi_chart_frame()`); it has no CDF, so a `cens()` or `trunc()`
#' response keeps the `\"xi\"` chart.
#'
#' @noRd
fam_gev_impl <- function(link, link_sigma, link_xi, chart) {
  lk <- mu_link(link, \"gen_extreme_value\")")
r("  frmtmb_family(
    \"gen_extreme_value\",
    accepts_aterms = c(\"weights\", \"cens\", \"trunc\"),
    dpars = c(\"mu\", \"sigma\", \"xi\"),
    links = list(mu = lk, sigma = lk_sigma, xi = lk_xi),
    lpdf = function(y, dpars, aterms) {
      ov <- osa_unwrap(y)
      g <- gev_parts(if (is.null(ov)) y else ov$y, dpars)",
"  core_lpdf <- function(y, dpars) {
    ov <- osa_unwrap(y)
    g <- gev_parts(if (is.null(ov)) y else ov$y, dpars)")
r("      ll <- -dpar_log(dpars, \"sigma\", lk_sigma) - g[[\"l\"]] - g[[\"h\"]] -
        exp_neg_capped(g[[\"h\"]]) - 1e8 * g[[\"out\"]] + xi_guard(dpars)
      if (is.null(ov)) ll else ll * ov$keep
    },
    lcdf = function(q, dpars, aterms) {
      exp(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]])) + xi_guard(dpars)
    },
    lccdf = function(q, dpars, aterms) {
      log(-expm1(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]]))) +
        xi_guard(dpars)
    },",
"    ll <- -dpar_log(dpars, \"sigma\", lk_sigma) - g[[\"l\"]] - g[[\"h\"]] -
      exp_neg_capped(g[[\"h\"]]) - 1e8 * g[[\"out\"]]
    if (is.null(ov)) ll else ll * ov$keep
  }
  lpdf <- function(y, dpars, aterms) core_lpdf(y, dpars) + xi_guard(dpars)
  lcdf <- function(q, dpars, aterms) {
    exp(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]])) + xi_guard(dpars)
  }
  lccdf <- function(q, dpars, aterms) {
    log(-expm1(-exp_neg_capped(gev_parts(q, dpars)[[\"h\"]]))) +
      xi_guard(dpars)
  }
  links <- list(mu = lk, sigma = lk_sigma, xi = lk_xi)
  if (identical(chart, \"tmp\")) {
    # no guard here: scale_xi() keeps every row inside the support, and
    # its range reaches below -1, as brms's does
    lpdf <- function(y, dpars, aterms) {
      frmtmb_gev_xi <- gev_scale_xi(dpars[[\"xi\"]], y, dpars[[\"mu\"]],
                                    dpars[[\"sigma\"]])
      RTMB::REPORT(frmtmb_gev_xi)
      dpars[[\"xi\"]] <- frmtmb_gev_xi
      core_lpdf(y, dpars)
    }
    lcdf <- NULL
    lccdf <- NULL
    links$xi <- get_link(\"identity\", dpar = \"xi\")
  }
  frmtmb_family(
    \"gen_extreme_value\",
    accepts_aterms = c(\"weights\", \"cens\", \"trunc\"),
    dpars = c(\"mu\", \"sigma\", \"xi\"),
    links = links,
    lpdf = lpdf,
    lcdf = lcdf,
    lccdf = lccdf,")
r("    post = list(
      fit_check = gev_fit_check,",
"    post = list(
      fit_check = gev_fit_check,
      # the sampler's chart for a scalar xi (gev_xi_chart_frame())
      tmp_xi_family = function() {
        fam_gev_impl(link, link_sigma, link_xi, chart = \"tmp\")
      },")
writeLines(strsplit(x, "\n")[[1]], p)
