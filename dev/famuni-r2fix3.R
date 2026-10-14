p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("gev_xi_chart_terms <- function(fit) {
  out <- list()
  for (rs in fit$spec$responses) {
    if (!identical(rs$family[[\"family\"]], \"gen_extreme_value\")) next
    r <- rs$resp_name
    key <- linpred_key(r, \"xi\")
    lp <- fit$frame[[\"linpreds\"]][[key]]
    if (is.null(lp) || !is.null(lp[[\"constant\"]]) ||
          !is.null(lp[[\"equate\"]]) || !is.null(lp[[\"nl_body\"]]) ||
          dpar_has_predictor(fit$spec, lp)) next
    av <- fit$frame[[\"aterm_values\"]][[r]] %||% list()
    blocked <- c(if (!is.null(av[[\"cens\"]])) \"cens()\",
                 if (!is.null(av[[\"trunc_lb\"]]) ||
                       !is.null(av[[\"trunc_ub\"]])) \"trunc()\")
    out[[length(out) + 1L]] <- list(resp = r, key = key, par = lp[[\"par\"]],
                                    idx = lp[[\"idx\"]], link = lp[[\"link\"]],
                                    blocked = blocked)
  }",
"gev_xi_chart_terms <- function(fit) {
  out <- list()
  for (rs in fit$spec$responses) {
    fname <- rs$family[[\"family\"]]
    # a mixture's GEV components carry xi1, xi2, ...; the chart covers a
    # plain GEV response only, and a mixture is listed so that the
    # sampler says so (a two-GEV mixture at xi -0.3 had 84 of 4000
    # transitions divergent, dev/famuni-rev2-out/mixsample.txt)
    mix <- !identical(fname, \"gen_extreme_value\")
    if (mix && !(startsWith(fname, \"mixture(\") &&
                   grepl(\"gen_extreme_value\", fname, fixed = TRUE))) next
    r <- rs$resp_name
    dps <- if (mix) grep(\"^xi[0-9]+$\", rs$family[[\"dpars\"]], value = TRUE)
           else \"xi\"
    for (dp in dps) {
      key <- linpred_key(r, dp)
      lp <- fit$frame[[\"linpreds\"]][[key]]
      if (is.null(lp) || !is.null(lp[[\"constant\"]]) ||
            !is.null(lp[[\"equate\"]]) || !is.null(lp[[\"nl_body\"]]) ||
            dpar_has_predictor(fit$spec, lp)) next
      av <- fit$frame[[\"aterm_values\"]][[r]] %||% list()
      blocked <- c(if (mix) \"a mixture()\",
                   if (!is.null(av[[\"cens\"]])) \"cens()\",
                   if (!is.null(av[[\"trunc_lb\"]]) ||
                         !is.null(av[[\"trunc_ub\"]])) \"trunc()\")
      out[[length(out) + 1L]] <- list(resp = r, dpar = dp, key = key,
                                      par = lp[[\"par\"]], idx = lp[[\"idx\"]],
                                      link = lp[[\"link\"]], blocked = blocked)
    }
  }")
r("#' `blocked`, the reasons the chart cannot be used (a `cens()` or
#' `trunc()` term, whose CDF the chart does not carry, or a second GEV
#' response, which would share the REPORTed name).",
"#' `blocked`, the reasons the chart cannot be used (a `cens()` or
#' `trunc()` term, whose CDF the chart does not carry, a second GEV
#' response, which would share the REPORTed name, or a [mixture()]
#' component, whose `xi1`, `xi2`, ... the chart does not map).")
# item 6: the stale comment
r("    # no guard here: scale_xi() keeps every row inside the support, and
    # its range reaches below -1, as brms's does",
"    # no guard here: gev_scale_xi() keeps every row inside the support
    # and holds its range above -1, where brms's goes on")
# item 1: the derivative claim
r("#' 0 where `x > 0` and `NaN` where `x <= 0`, with derivatives 0 where it
#' is 0, on the RTMB tape as on doubles: a term that takes a parameter
#' out of the likelihood where the density does not exist.",
"#' 0 where `x > 0` and `NaN` where `x <= 0`, on the RTMB tape as on
#' doubles: a term that takes a parameter out of the likelihood where
#' the density does not exist. Its derivatives are 0 where it is 0, with
#' three exceptions at the ends of the doubles (the review's round 2,
#' dev/famuni-rev2-out/guard.txt and hess.txt): the gradient is NaN for
#' `0 < x < 2.3e-308`, where `x` is subnormal; the Hessian is NaN for
#' `0 < x <= 1e-155`; and the value itself is NaN for `x >= 9e307`, where
#' `u` overflows. The GEV and the Frechet cannot reach these: `1 + xi`
#' and `nu - 1` are 0 or at least 1.1e-16 in double precision. frmtmb.eam's
#' wiener guard reaches them only with a parameter below 1e-154.")
writeLines(strsplit(x, "\n")[[1]], p)
