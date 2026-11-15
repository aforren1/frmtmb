# Refusal of an asymmetric Laplace model whose Laplace approximation has
# no maximum in the scale of its random effects.

# The families whose log density is piecewise linear in mu: a kink at
# mu and a straight line on each side, so its second derivative in mu
# is 0 wherever it is defined. A zero of zero_inflated_asym_laplace
# reads only zi.
ald_laplace_families <- c("asym_laplace", "zero_inflated_asym_laplace")

# The built-in structures whose block density is gaussian with a
# covariance that a free scale multiplies (their sd, sds or sdgp, all
# of them for a block with several), and whose Laplace objective was
# seen to rise along that scale (dev/wcheck2-ald-profile.R). The other
# scaled structures (hetar1, ou, toep, homtoep, homcs, exp, gau, mat,
# rr, gr_prec) were not profiled along their scale, and the exact
# gp()'s objective fell at large sd where TMB's inner search failed, so
# they stay open. equalto() fixes the covariance, and a Student-t block
# is not gaussian.
ald_scaled_covstructs <- c("us", "diag", "homdiag", "cs", "ar1", "gr_cov",
                           "smooth", "hsgp")

# The prior classes that leave a block's scale free: a prior on any
# other class (sd, sds, sdgp, or a structure's own parameters) can bound
# the scale, so the check stays open there
ald_scale_free_classes <- c("b", "Intercept", "sigma", "quantile", "zi",
                            "cor", "L", "lscale")

#' Refuse an asym_laplace() or zero_inflated_asym_laplace() model whose
#' random effects in mu have no estimate under the Laplace approximation.
#'
#' With mu linear in the random effects b (the identity link, no
#' nonlinear body), the log density is piecewise linear in b, so the
#' data add no curvature: the Hessian of the inner problem is the
#' prior's alone, `Sigma^-1`. The approximation is then
#' `max_b (l(b) + log N(b; 0, Sigma))` plus
#' `q/2 log(2 pi) - 1/2 log det H`, which is
#' `max_b (l(b) - b' Sigma^-1 b / 2)`: the normalizing constant of the
#' prior and the curvature term cancel. Multiplying a block's covariance
#' by c^2 > 1 multiplies its penalty by 1 / c^2, so the approximation
#' never falls as the
#' block's scale grows: it has no maximum at a finite scale (or, where
#' it stays flat, no unique one). It holds for every value of the other
#' parameters, so a prior on them changes nothing; a prior on the scale
#' itself can bound it, and the check then stays open.
#'
#' The proof needs the zero curvature, so the check stays open for:
#' another link (`rho'(u) f''(eta)` enters the curvature), a nonlinear
#' body, `cens()` or `trunc()` (a tail probability is curved), a random
#' effect of the same response in sigma, quantile or zi (the cross
#' derivatives with mu are not 0), a block shared with another response,
#' `mi()` or `me()` latent values, a residual correlation term, a
#' Student-t block or a structure another package registered anywhere
#' in that mu (its curvature depends on its mode, or is not known, so no
#' block of that mu is refused), a fixed-covariance block (not refused
#' itself), and `quadrature = TRUE` or `importance` (the marginal is
#' then an integral, smooth in the scale). Measured on rellib-r21
#' (dev/reviews/wcheck2-review-al.R): the approximation at exact inner
#' modes rose along log sd in all 30 profiles, and TMB's own inner
#' search left modes up to 2.06 away from them.
#'
#' @noRd
check_ald_laplace <- function(spec, frame, prior, exact = FALSE) {
  if (isTRUE(exact) || length(frame[["mi_map"]] %||% list()) ||
        length(frame[["me"]] %||% list()) ||
        length(frame[["autocor"]] %||% list())) {
    return(invisible())
  }
  if (length(prior)) {
    if (!inherits(prior, "frmtmb_priorlist")) return(invisible())
    cls <- vapply(prior, function(s) s$class %||% "", "")
    if (!all(cls %in% ald_scale_free_classes)) return(invisible())
  }
  for (r in names(spec$responses)) {
    hit <- ald_laplace_blocks(spec, frame, r)
    if (length(hit)) {
      frm_stop(ald_laplace_message(spec, r, hit), call. = FALSE,
               class = "frmtmb_ald_no_laplace")
    }
  }
  invisible()
}

#' The blocks of response `r` whose scale the Laplace approximation
#' cannot estimate (see `check_ald_laplace()`), as their terms; empty
#' when the response is not decided.
#'
#' @noRd
ald_laplace_blocks <- function(spec, frame, r) {
  resp <- spec$responses[[r]]
  fam <- resp$family[["family"]] %||% ""
  if (!fam %in% ald_laplace_families) return(character())
  av <- frame[["aterm_values"]][[r]]
  if (any((av[["cens"]] %||% 0) != 0) || !is.null(av[["trunc_lb"]]) ||
        !is.null(av[["trunc_ub"]])) {
    return(character())
  }
  key <- linpred_key(r, "mu")
  lp <- frame[["linpreds"]][[key]]
  if (is.null(lp) || !is.null(lp[["nl_body"]]) ||
        !identical(lp[["link"]][["name"]] %||% "", "identity")) {
    return(character())
  }
  # the response's blocks: one with a component off its mu couples mu's
  # effects with a curved density, and the proof does not hold
  keys <- lapply(frame[["re_blocks"]], function(bk) {
    vapply(bk[["components"]] %||% list(), function(cp) {
      cp[["lp_key"]] %||% ""
    }, "")
  })
  dpars <- names(resp$family[["links"]] %||% list())
  own <- vapply(keys, function(k) {
    any(k %in% vapply(dpars, function(dp) linpred_key(r, dp), ""))
  }, NA)
  if (!any(own)) return(character())
  in_mu <- vapply(keys, function(k) length(k) && all(k == key), NA)
  if (any(own & !in_mu)) return(character())
  # the proof needs H = Sigma^-1 for the whole inner problem: a
  # Student-t block adds a curvature that depends on its mode, which
  # moves with every other block's sd, and another package's block
  # states no density here, so beside either nothing is decided
  curved <- vapply(frame[["re_blocks"]], function(bk) {
    !identical(bk[["covstruct_package"]] %||% "frmtmb", "frmtmb") ||
      !is.null(bk[["dist_nu"]]) || !is.null(bk[["t_share"]]) ||
      bk[["covstruct"]] %in% c("us_t", "diag_t")
  }, NA)
  if (any(own & curved)) return(character())
  multi <- length(spec$responses) > 1L
  terms <- character()
  for (i in which(own)) {
    bk <- frame[["re_blocks"]][[i]]
    if (!identical(bk[["covstruct_package"]] %||% "frmtmb", "frmtmb") ||
          !(bk[["covstruct"]] %in% ald_scaled_covstructs) ||
          !is.null(bk[["dist_nu"]]) || !length(bk[["theta_idx"]])) {
      next
    }
    lab <- bk[["term_label"]] %||% ""
    pre <- paste0(resp$resp_name %||% r, " ")
    if (multi && startsWith(lab, pre)) lab <- substring(lab, nchar(pre) + 1L)
    terms <- c(terms, if (bk[["covstruct"]] %in% c("smooth", "hsgp")) {
      # a smooth's or a gp's label is its term already: "s(z)", "gp(z)"
      lab
    } else if (identical(bk[["covstruct"]], "gr_cov")) {
      paste0("(", lab, ") with gr(cov = )")
    } else if (multi) {
      block_term_in_response(bk, r, resp$resp_name %||% r)
    } else block_term_as_written(bk))
  }
  unique(terms)
}

#' The refusal's text.
#'
#' @noRd
ald_laplace_message <- function(spec, r, terms) {
  fam <- spec$responses[[r]]$family[["family"]]
  multi <- length(spec$responses) > 1L
  paste0("frm() cannot fit ", paste0("`", terms, "`", collapse = ", "),
         if (multi) paste0(" in the mu of response ", r), " under the ",
         "Laplace approximation: the ", fam, " log density is linear in ",
         "mu on each side of its mode, so the data give these random ",
         "effects no curvature, and the approximation never falls as ",
         "their standard deviation grows, so it has no maximum. ",
         "frmtmb.sample::frm_sample() samples the model without the ",
         "approximation. quadrature = TRUE integrates one grouping ",
         "factor exactly, but on 40 test data sets its rule broke down on ",
         "24.")
}
