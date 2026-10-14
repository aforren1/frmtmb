p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
anchor <- "# ------------------------------------------------------- discrete weibull"
stopifnot(grepl(anchor, x, fixed = TRUE))
add <- r"--(#' The largest or smallest element of `v`, on the tape: a pairwise
#' tree of the tape-safe `ad_max()` or `ad_min()`, `log2(n)` vectorized
#' steps rather than `n` scalar ones. An odd length repeats its last
#' element, which changes neither.
#'
#' @noRd
ad_tree_extreme <- function(v, f) {
  n <- length(v)
  while (n > 1L) {
    if (n %% 2L) v <- v[c(seq_len(n), n)]
    i <- seq.int(1L, length(v), by = 2L)
    v <- f(v[i], v[i + 1L])
    n <- length(v)
  }
  v
}

#' brms's `scale_xi()`: the GEV shape from an unbounded `tmp_xi`,
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
}

#' The sampler's chart for a scalar GEV shape
#'
#' A `gen_extreme_value()` `xi` without a formula is sampled as brms
#' samples it, as `tmp_xi` on the whole line, mapped into the support the
#' current `mu` and `sigma` leave by brms's `scale_xi()`. Sampled on its
#' `log1p` scale, the support's edge is a wall in the sampled space: on
#' the review's data (n = 40, `y ~ 1`, 4 chains x 2000) 3388 of 4000
#' transitions diverged at xi = -0.85, and none did through `tmp_xi`
#' (dev/famuni-out/tmpxi.txt).
#'
#' `gev_xi_chart_terms()` lists each such `xi`: its response, linear
#' predictor key, parameter component and index, the fit's link, and
#' `blocked`, the reasons the chart cannot be used (a `cens()` or
#' `trunc()` term, whose CDF the chart does not carry, or a second GEV
#' response, which would share the REPORTed name). `gev_xi_chart_frame()`
#' returns the frame with each unblocked one in the `tmp_xi` chart, and
#' `gev_xi_chart_start()` carries a parameter list's `xi` coefficients
#' to `tmp_xi` at the fit's `mu` and `sigma`. `dpar_has_predictor()`
#' tells `xi ~ 1`, which brms samples on `log1p`, from no formula.
#'
#' @param fit A `frmtmb_fit`.
#' @param frame A model frame.
#' @param terms What `gev_xi_chart_terms()` returned.
#' @param pars A parameter list in the frame's template.
#' @param spec,lp A model spec and one of its frame's linear predictors.
#' @keywords internal
#' @name frmtmb-gev-xi-chart
#' @rawNamespace export(gev_xi_chart_terms, gev_xi_chart_frame,
#'   gev_xi_chart_start, dpar_has_predictor)
NULL

#' @rdname frmtmb-gev-xi-chart
gev_xi_chart_terms <- function(fit) {
  out <- list()
  for (rs in fit$spec$responses) {
    if (!identical(rs$family[["family"]], "gen_extreme_value")) next
    r <- rs$resp_name
    key <- linpred_key(r, "xi")
    lp <- fit$frame[["linpreds"]][[key]]
    if (is.null(lp) || !is.null(lp[["constant"]]) ||
          !is.null(lp[["equate"]]) || !is.null(lp[["nl_body"]]) ||
          dpar_has_predictor(fit$spec, lp)) next
    av <- fit$frame[["aterm_values"]][[r]] %||% list()
    blocked <- c(if (!is.null(av[["cens"]])) "cens()",
                 if (!is.null(av[["trunc_lb"]]) ||
                       !is.null(av[["trunc_ub"]])) "trunc()")
    out[[length(out) + 1L]] <- list(resp = r, key = key, par = lp[["par"]],
                                    idx = lp[["idx"]], link = lp[["link"]],
                                    blocked = blocked)
  }
  if (length(out) > 1L) {
    out <- lapply(out, function(tm) {
      tm$blocked <- c(tm$blocked, "a second gen_extreme_value() response")
      tm
    })
  }
  out
}

#' @rdname frmtmb-gev-xi-chart
gev_xi_chart_frame <- function(frame, terms) {
  for (tm in terms) {
    if (length(tm$blocked)) next
    fam <- frame[["spec"]]$responses[[tm$resp]]$family
    frame[["spec"]]$responses[[tm$resp]]$family <- fam[["post"]][["tmp_xi_family"]]()
    frame[["linpreds"]][[tm$key]][["link"]] <- get_link("identity", dpar = "xi")
  }
  frame
}

#' @rdname frmtmb-gev-xi-chart
gev_xi_chart_start <- function(fit, terms, pars) {
  multi <- length(fit$spec$responses) > 1L
  for (tm in terms) {
    if (length(tm$blocked)) next
    rs <- if (multi) tm$resp
    lpr <- function(dp) frm_linpred(fit, dpar = dp, type = "response", resp = rs)
    x <- (fit$frame[["y"]][[tm$resp]] - lpr("mu")) / lpr("sigma")
    b <- -1 / range(x)
    xi <- lpr("xi")[1L]
    q <- (xi - min(b)) / (max(b) - min(b))
    # the maximum is inside the support, so q is inside (0, 1); held off
    # its ends so a start on the edge has a finite tmp_xi
    pars[[tm$par]][tm$idx] <- stats::qlogis(min(max(q, 1e-8), 1 - 1e-8))
  }
  pars
}

)--"
x <- sub(anchor, paste0(add, anchor), x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
