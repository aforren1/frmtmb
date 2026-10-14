p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new, all = FALSE) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- if (all) gsub(old, new, x, fixed = TRUE) else
    sub(old, new, x, fixed = TRUE)
}
r("mo_sample_frame <- function(fit) {
  if (!length(simplex_frame_terms(fit))) return(fit$frame)
  mo_chart_frame(fit$frame, \"softmax\")
}",
"mo_sample_frame <- function(fit) {
  fr <- if (!length(simplex_frame_terms(fit))) fit$frame else {
    mo_chart_frame(fit$frame, \"softmax\")
  }
  # a scalar GEV xi in brms's tmp_xi chart (frmtmb's gev_xi_chart_frame())
  gev_xi_chart_frame(fr, gev_xi_chart_terms(fit))
}")
r("mo_sample_pars <- function(fit, pars) {
  for (tm in simplex_frame_terms(fit)) {
    pars[[tm$zeta]] <- mo_coords(mo_simplex(pars[[tm$zeta]]), \"softmax\")
  }
  pars
}",
"mo_sample_pars <- function(fit, pars) {
  for (tm in simplex_frame_terms(fit)) {
    pars[[tm$zeta]] <- mo_coords(mo_simplex(pars[[tm$zeta]]), \"softmax\")
  }
  gev_xi_chart_start(fit, gev_xi_chart_terms(fit), pars)
}")
r("mo_sample_entries <- function(entries) {
  lapply(entries, function(e) {
    if (identical(e$dist$kind, \"dirichlet\")) e$dist$chart <- \"softmax\"
    e
  })
}",
"mo_sample_entries <- function(entries, fit) {
  terms <- Filter(function(tm) !length(tm$blocked), gev_xi_chart_terms(fit))
  lapply(entries, function(e) {
    if (identical(e$dist$kind, \"dirichlet\")) e$dist$chart <- \"softmax\"
    # a prior on a scalar xi is brms's prior on tmp_xi, the coefficient
    # the chart samples, so it is read on that coefficient as it stands
    for (tm in terms) {
      if (identical(e$comp, tm$par) && length(e$idx) &&
            all(e$idx %in% tm$idx)) {
        e$scale <- \"tmp_xi\"
        e$link <- NULL
      }
    }
    e
  })
}")
r("  entries <- mo_sample_entries(entries)", "  entries <- mo_sample_entries(entries, fit)", all = TRUE)
# the route: a scalar-xi GEV is retaped in its chart like a mo() simplex
r("  } else if (length(simplex_frame_terms(fit))) {",
  "  } else if (length(simplex_frame_terms(fit)) || gev_chart) {")
r("  ncp <- ncp_plan(fit, reparameterize, laplace, ri$entries)",
  "  # a scalar GEV xi is sampled through brms's tmp_xi where it can be,
  # and said so where it cannot (gev_xi_blocked_text())
  gev_terms <- gev_xi_chart_terms(fit)
  gev_chart <- any(vapply(gev_terms, function(tm) !length(tm$blocked), NA))
  for (tm in gev_terms) {
    if (length(tm$blocked)) {
      frm_warning(\"frm_sample(): \", gev_xi_blocked_text(tm, \"\"),
                  call. = FALSE)
    }
  }
  ncp <- ncp_plan(fit, reparameterize, laplace, ri$entries)")
r("  if (mo_soft) m <- mo_draws_to_chart(m, fit, laplace)",
  "  if (mo_soft) m <- mo_draws_to_chart(m, fit, laplace)
  if (gev_chart && !identical(obj, fit$obj)) {
    m <- gev_draws_to_xi(m, obj, fit, gev_terms, laplace)
  }")
# helpers
r("#' The draws matrix with each `mo()` simplex's softmax columns carried to",
"#' Why a scalar GEV `xi` is not sampled through `tmp_xi`, and what that
#' costs, for the warning and the disclosure note.
#'
#' @noRd
gev_xi_blocked_text <- function(tm, of = \"\") {
  paste0(\"xi\", of, \" of gen_extreme_value() is sampled on its \",
         tm$link$name, \" scale, not through brms's tmp_xi, because of \",
         paste(tm$blocked, collapse = \" and \"), \". The support's edge \",
         \"is then a wall the sampler meets: on the data of \",
         \"dev/famuni-out/tmpxi.txt, 957 to 3388 of 4000 transitions \",
         \"diverged that way and none through tmp_xi. Check the \",
         \"divergences, or give xi a formula (xi ~ 1)\")
}

#' The draws matrix with a scalar GEV `xi` carried from `tmp_xi`, the
#' chart the sampler ran in, to the fit's own coefficient: each draw's
#' `xi` is the one its density used, which the tape REPORTs, under the
#' fit's link. brms's `compute_xi()` does the same after sampling.
#'
#' @noRd
gev_draws_to_xi <- function(m, obj, fit, terms, laplace) {
  terms <- Filter(function(tm) !length(tm$blocked), terms)
  if (!length(terms)) return(m)
  tm <- terms[[1L]]
  pidx <- draws_par_index(fit, laplace = laplace)
  tpl <- fit$frame[[\"par_template\"]]
  fx <- if (identical(tm$par, \"betad\")) fit$frame[[\"betad_fixed_idx\"]]
  free <- setdiff(seq_along(tpl[[tm$par]]), fx)
  j <- pidx[[tm$par]][match(tm$idx, free)]
  np <- if (laplace) length(obj$par) else length(obj$env$last.par)
  for (i in seq_len(nrow(m))) {
    row <- as.numeric(m[i, seq_len(np)])
    rep <- if (laplace) {
      obj$fn(row)
      obj$report()
    } else {
      obj$report(row)
    }
    m[i, j] <- tm$link$linkfun(rep[[\"frmtmb_gev_xi\"]][1L])
  }
  m
}

#' The draws matrix with each `mo()` simplex's softmax columns carried to")
writeLines(strsplit(x, "\n")[[1]], p)
