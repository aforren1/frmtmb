# Release fixes of the 0.77.0 consolidation, sourced by
# dev/rel077-resolve.R (local = TRUE, so rd(), wr(), edit_file() and
# restore() are the resolver's), after the textual merge and the hand
# resolutions. Each edit checks that its old text is there exactly once,
# so a merge rerun applies them all or stops at the first that no longer
# fits. They come from the lanes' reviews (dev/reviews/2026-10-27-<lane>.md,
# the short checks of the fix rounds), named by the coordinator as fixes
# for the release. NEWS edits of fact are in dev/rel077-news-split.R,
# because the NEWS merge reads the lanes' own text after this file runs.
# The test blocks this file appends are dev/rel077-fix-*-test.R.
#
#   (sourced; not run on its own)
append_test <- function(f, blk, marker) {
  x <- rd(f)
  stopifnot(!any(grepl(marker, x, fixed = TRUE)))
  wr(c(x, "", rd(blk)), f)
  cat("appended", blk, "to", f, "\n")
}

# ---- S1 (sigmainn). insight::get_sigma() keeps insight's default --------
# sigmainn's review, short check of fix round 1, S1: the lane's
# .get_sigma.frmtmb_fit returned sigma() itself, so on a plain
# bf(y ~ x, sigma ~ x2) fit, where sigma() is NA, insight's fallback
# sqrt(deviance / df) (1.7285282 on rellib-r14, seed 3) became NA. The
# method only has to refuse a fit saved with the old sigma, then hand
# over (dev/review-sigmainn-nextmethod.R in the lane's tree). The NEWS
# line is reworded in dev/rel077-news-split.R.
edit_file("R/insight.R",
"# insight::get_sigma() is not generic; it dispatches on its internal
# .get_sigma(), whose default falls back on sqrt(deviance / df) when
# sigma() stops. For a fit saved with an old sigma under ar(cov = TRUE)
# that fallback was 1.23, neither sd. This method gives sigma() itself,
# as the default does when sigma() answers, and lets sigma()'s refusal
# reach the caller.
#' @exportS3Method insight::.get_sigma
.get_sigma.frmtmb_fit <- function(x, verbose = TRUE, ...) {
  s <- stats::sigma(x)
  class(s) <- c(\"insight_aux\", class(s))
  s
}",
"# insight::get_sigma() is not generic; it dispatches on its internal
# .get_sigma(), whose default reads sigma() and falls back on
# sqrt(deviance / df) when sigma() gives nothing. For a fit saved with an
# old sigma under ar(cov = TRUE) that fallback was 1.23, neither sd, so
# such a fit is refused here; any other fit gets insight's default, NA
# from sigma() under sigma ~ x included.
#' @exportS3Method insight::.get_sigma
.get_sigma.frmtmb_fit <- function(x, verbose = TRUE, ...) {
  autocor_check_fit(x)
  NextMethod()
}")
append_test("tests/testthat/test-autocor-sigma.R",
            "dev/rel077-fix-getsigma-test.R", "insight's default on a")

# ---- R1 (reslot). A fit-route table through as.data.frame() -------------
# reslot's review, short check of fix round 2, R1: as.data.frame() drops
# a default_prior() table's route (by design: the coerced table is the
# plain data frame a reader would write, test-get-prior-route.R), so
# frm_sample(prior = as.brmsprior(as.data.frame(pe))) read the fit-route
# (flat) rows as brms's (flat) and sampled flat on the intercept, sigma,
# sd and cor in silence, the outcome the fit-route refusal exists to
# prevent. Keeping the mark through as.data.frame() is not clean, so:
# as.brmsprior() of a plain data frame that holds frmtmb's printed
# "(flat)" on a row whose source is "default" (brms's own tables store
# "" there) marks the priorlist "table", the mark travels as "fit" does,
# and frm_sample() refuses such a prior where its (flat) rows replace a
# proper brms default, which a route = "sample" table never does. The
# plain-data-frame refusal says to pass the table as it is.
pf <- "R/priors.R"
edit_file(pf,
"#' A priorlist marked with the route of the default_prior() table it was
#' read from. A route = \"fit\" table's (flat) rows are frm()'s defaults,
#' not a user's choice, and frmtmb.sample refuses it rather than read
#' them as brms's (flat) (`sample_resolve_priors()`); the mark travels
#' through `+`, `c()` and a column edit.
#'
#' @noRd
with_prior_route <- function(pl, route) {
  if (!is.null(pl) && identical(route, \"fit\")) attr(pl, \"route\") <- route
  pl
}

#' The route mark of several priorlists combined: \"fit\" when any of them
#' carries it.
#'
#' @noRd
combined_prior_route <- function(parts) {
  r <- vapply(parts, function(p) attr(p, \"route\") %||% \"\", \"\")
  if (any(r == \"fit\")) \"fit\"
}",
"#' A priorlist marked with the route of the default_prior() table it was
#' read from. A route = \"fit\" table's (flat) rows are frm()'s defaults,
#' not a user's choice, and frmtmb.sample refuses it rather than read
#' them as brms's (flat) (`sample_resolve_priors()`); the mark travels
#' through `+`, `c()` and a column edit. The mark \"table\" is a plain
#' data frame of a frmtmb table, whose route as.data.frame() dropped:
#' frmtmb.sample refuses it where its (flat) rows replace a proper brms
#' default, which only a route = \"fit\" table's rows do.
#'
#' @noRd
with_prior_route <- function(pl, route) {
  if (!is.null(pl) && length(route) == 1L && route %in% c(\"fit\", \"table\")) {
    attr(pl, \"route\") <- route
  }
  pl
}

#' The route mark of several priorlists combined: \"fit\" when any of them
#' carries it, else \"table\" when any carries that.
#'
#' @noRd
combined_prior_route <- function(parts) {
  r <- vapply(parts, function(p) attr(p, \"route\") %||% \"\", \"\")
  if (any(r == \"fit\")) \"fit\" else if (any(r == \"table\")) \"table\"
}")
edit_file(pf,
"  route <- if (inherits(x, \"frmtmb_prior_rows\")) attr(x, \"route\")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!\"prior\" %in% names(x)) {
    frm_stop(\"as.brmsprior() needs a `prior` column\", call. = FALSE)
  }",
"  route <- if (inherits(x, \"frmtmb_prior_rows\")) attr(x, \"route\")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!\"prior\" %in% names(x)) {
    frm_stop(\"as.brmsprior() needs a `prior` column\", call. = FALSE)
  }
  # frmtmb's own table, coerced: its route is gone, and its (flat) rows
  # may be frm()'s defaults; brms's tables store \"\" there, not \"(flat)\"
  if (is.null(route) && \"source\" %in% names(x) &&
      any(x[[\"prior\"]] %in% \"(flat)\" & x[[\"source\"]] %in% \"default\")) {
    route <- \"table\"
  }")
edit_file(pf,
"      \"prior = got a plain data frame. brms refuses one too (\\\"Argument \",
      \"'prior' must be a 'brmsprior' object\\\"); turn it into a prior \",
      \"with as.brmsprior(), which reads brms's columns (prior, class, \",
      \"coef, group, resp, dpar, nlpar, lb, ub)\"))",
"      \"prior = got a plain data frame. brms refuses one too (\\\"Argument \",
      \"'prior' must be a 'brmsprior' object\\\"); turn it into a prior \",
      \"with as.brmsprior(), which reads brms's columns (prior, class, \",
      \"coef, group, resp, dpar, nlpar, lb, ub). A default_prior() table \",
      \"is passed as it is, since as.data.frame() drops its route; for \",
      \"frmtmb.sample::frm_sample(), take it with route = \\\"sample\\\"\"))")
sa <- "extensions/frmtmb.sample/R/sample.R"
edit_file(sa,
"  if (!is.null(attr(prior, \"route\"))) attr(prior, \"route\") <- NULL
  st <- prior_stack(base, prior)",
"  # a plain data frame of a frmtmb table lost its route: checked below,
  # against what its (flat) rows replace
  table_mark <- identical(attr(prior, \"route\"), \"table\")
  if (!is.null(attr(prior, \"route\"))) attr(prior, \"route\") <- NULL
  st <- prior_stack(base, prior)")
edit_file(sa,
"  defs <- if (defaults) default_priors_for(fit)
  if (!is.null(defs) && !is.null(user_pl)) {",
"  defs <- if (defaults) default_priors_for(fit)
  defs_all <- defs
  if (!is.null(defs) && !is.null(user_pl)) {")
edit_file(sa,
"  ri <- if (!is.null(eff)) resolve_prior_input(fit, eff)
  # announced even when there are no defaults at all",
"  ri <- if (!is.null(eff)) resolve_prior_input(fit, eff)
  if (table_mark && defaults) {
    lost <- table_flat_lost(fit, ri, defs_all)
    if (length(lost)) {
      frm_stop(\"prior = was read from a default_prior() or validate_prior() \",
               \"table that as.data.frame() turned into a plain data frame, \",
               \"which drops its route, and its (flat) rows would replace \",
               \"brms's default on \", paste(lost, collapse = \", \"),
               \" with (flat). A table of route = \\\"fit\\\", the default, \",
               \"lists frm()'s defaults as (flat). For brms's defaults with \",
               \"your rows edited in, pass the table itself or start from \",
               \"default_prior(..., route = \\\"sample\\\"); for brms's (flat) \",
               \"in every slot, write prior = \\\"flat\\\"; or pass only the \",
               \"specifications you set\", call. = FALSE)
    }
  }
  # announced even when there are no defaults at all")
edit_file(sa,
"#' What a user's specification of a slot keeps from brms's default
#' there, since brms keeps both.",
"#' The classes of brms's proper defaults whose every position the
#' effective prior `ri` leaves (flat): what the (flat) rows of a coerced
#' frmtmb table replaced. A default that is already flat (on an interval,
#' a flat LKJ) or dirichlet(1), which brms's (flat) on a simplex is, is
#' no loss.
#'
#' @noRd
table_flat_lost <- function(fit, ri, defs) {
  is_flat <- function(e) {
    identical(e[[\"dist\"]][[\"kind\"]], \"flat\") ||
      isTRUE(e[[\"dist\"]][[\"flat\"]])
  }
  key <- function(es) {
    unlist(lapply(es, function(e) paste(e[[\"comp\"]], e[[\"idx\"]])))
  }
  flat <- key(Filter(is_flat, ri[[\"entries\"]] %||% list()))
  out <- character(0)
  for (d in unclass(defs %||% list())) {
    dd <- d[[\"dist\"]]
    if (is.null(dd) || isTRUE(d[[\"flat_interval\"]]) ||
          dd[[\"kind\"]] %in% c(\"flat\", \"dirichlet\") ||
          isTRUE(dd[[\"flat\"]])) next
    di <- resolve_prior_input(fit, structure(list(d),
                                             class = \"frmtmb_priorlist\"))
    k <- key(di[[\"entries\"]])
    if (length(k) && all(k %in% flat)) out <- c(out, spec_spelling(d)$class)
  }
  unique(out)
}

#' What a user's specification of a slot keeps from brms's default
#' there, since brms keeps both.")
append_test("extensions/frmtmb.sample/tests/testthat/test-flat-prior.R",
            "dev/rel077-fix-route-test.R", "as.data.frame() drops the route")
# reslot's own assertion pinned the trap as \"read as written\"; a table
# built by hand (no `source` column) still is
edit_file("extensions/frmtmb.sample/tests/testthat/test-flat-prior.R",
"  # a table with no route (built by hand) is read as written: its
  # (flat) rows are brms's (flat)
  hand <- as.brmsprior(as.data.frame(p))
  expect_null(attr(hand, \"route\"))
  eh <- fp_rp(fit, hand)$ri$entries
  expect_true(\"theta 1 flat\" %in% key(eh))",
"  # a table with no route (built by hand, so with no `source` column)
  # is read as written: its (flat) rows are brms's (flat)
  hand <- as.brmsprior(as.data.frame(p)[, setdiff(names(p), \"source\")])
  expect_null(attr(hand, \"route\"))
  eh <- fp_rp(fit, hand)$ri$entries
  expect_true(\"theta 1 flat\" %in% key(eh))")

# ---- R2 (reslot). A list of specifications ------------------------------
# reslot's review, short check of fix round 2, R2: frm(prior = list(sp,
# sp2)) got the generic "not a list of length 2"; brms combines with +
edit_file(pf,
"  paste0(\"`prior` must be a set_prior() specification, a brms prior or a \",
         \"default_prior() table, not \", arg_desc(x))",
"  if (is.list(x) && length(x) && !is.data.frame(x) &&
        all(vapply(x, inherits, NA, \"frmtmb_priorlist\"))) {
    return(paste0(
      \"prior = got a list of \", length(x), \" set_prior() specifications, \",
      \"which frmtmb does not read as one prior. Combine them with + or \",
      \"c(), as in brms: sp1 + sp2, or c(sp1, sp2)\"))
  }
  paste0(\"`prior` must be a set_prior() specification, a brms prior or a \",
         \"default_prior() table, not \", arg_desc(x))")
append_test("tests/testthat/test-prior-one-slot.R",
            "dev/rel077-fix-priorlist-test.R", "a list of specifications")

# ---- S1, S2 (mointerop). ?frmtmb-marginaleffects ----------------------
# mointerop's review, short check of fix round 1: HC1's k counts the
# determined steps of each mo() simplex too, and drops a one-value term
# (interop_hc1_k(), which matches lm()); frmtmb does not compute
# leverages, which are defined for a plain gaussian or glm fit
mi <- "R/interop.R"
edit_file(mi,
"#'     `lm()`, with `k` the number of population-level coefficients
#'     (those of every linear predictor but the distributional
#'     parameters'). marginaleffects maps `\"stata\"` to HC1 for `lm()`.}",
"#'     `lm()`, with `k` the number of population-level coefficients
#'     (those of every linear predictor but the distributional
#'     parameters') plus the steps of each `mo()` simplex that the data
#'     determine; a `mo()` term on one value counts nothing, as `lm()`
#'     does not count a column it aliases. So `k` is the column count of
#'     `lm()` on `factor()` of the `mo()` variable. marginaleffects maps
#'     `\"stata\"` to HC1 for `lm()`.}")
edit_file(mi,
"#' marginaleffects maps to HC3) weight each row by its leverage, which a
#' fit by maximum likelihood with random effects or distributional
#' parameters does not define.",
"#' marginaleffects maps to HC3) weight each row by its leverage, the
#' diagonal of a hat matrix, which frmtmb does not compute.")
edit_file(mi,
"    paste0(\"It weights each row by its leverage, the diagonal of a hat \",
           \"matrix, which a model fitted by maximum likelihood with a \",
           \"Laplace approximation (random effects, distributional \",
           \"parameters) does not define\",",
"    paste0(\"It weights each row by its leverage, the diagonal of a hat \",
           \"matrix, which frmtmb does not compute\",")
edit_file(mi,
"#'   `n / (n - k)`, with `k` the population-level coefficients
#'   (`interop_hc1_k()`).",
"#'   `n / (n - k)`, with `k` the population-level coefficients and the
#'   determined steps of each `mo()` simplex (`interop_hc1_k()`).")

# ---- priordraws. The tests its review asks for on the merge -------------
# priordraws' review, finding 1 and section 13 of its findings: on the
# merge with reslot, default ar/ma (cov = TRUE) and cosy get no prior
# column and "only" is refused naming them; finding 2: under sigmainn,
# "only" sigma under ar(cov = TRUE) is its prior; mointerop: a two-value
# mo() under "yes" draws its one weight as 1; reslot: under
# prior = "flat" the flat entries get no prior draws, the evidence ratio
# is NA, and "only" is refused. dev/rel077-overlap.R is the probe.
pdt <- "dev/rel077-fix-priordraws-test.R"
if (file.exists(pdt)) {
  append_test("extensions/frmtmb.sample/tests/testthat/test-prior-draws.R",
              pdt, "the 0.77.0 consolidation")
} else cat("NOTE: no", pdt, "yet\n")

# ---- C1 (consolidation). brms's (flat) on a mo() or family simplex ----
# Found by the mointerop x reslot probe (dev/rel077-overlap.R mr,
# dev/rel077-mr-probe.R): frm_sample(prior = set_prior("", class =
# "simo", coef = "mok1")) stopped with "Invalid Dirichlet prior. Expected
# input of length 3 ... got 0", on reslot's lane too, with any number of
# categories. A specification with no density has no concentration to
# check; the resolver then skips it (reslot), so the simplex keeps its
# default dirichlet(1), which is brms's reading of "" there: brms keeps
# dirichlet_lpdf(simo_1 | con_simo_1) with con_simo_1 = 1 under "".
edit_file(pf,
"      a <- s$dist$alpha
      if (length(a) != 1L && length(a) != mi$D) {",
"      a <- s$dist$alpha
      if (!is.null(s$dist) && length(a) != 1L && length(a) != mi$D) {")
edit_file(pf,
"    a <- s$dist$alpha
    if (length(a) != 1L && length(a) != tm$D) {",
"    a <- s$dist$alpha
    if (!is.null(s$dist) && length(a) != 1L && length(a) != tm$D) {")
append_test("extensions/frmtmb.sample/tests/testthat/test-flat-prior.R",
            "dev/rel077-fix-flat-overlap-test.R",
            "(flat) meets sigmainn and mointerop")

# ---- the mointerop x sigmainn overlap, as a test ------------------------
append_test("tests/testthat/test-autocor-sigma.R",
            "dev/rel077-fix-me-sigma-test.R",
            "marginaleffects reads sigma as fitted() does")

# ---- after the release review: F1 to F6 ---------------------------------
# dev/reviews/2026-10-28-release.md (FIX, one major); the fixes and their
# reasons are dev/rel077-release-fix2.R
source("dev/rel077-release-fix2.R", local = TRUE)
