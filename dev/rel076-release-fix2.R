# Release fixes after the release review (dev/reviews/2026-10-27-release.md,
# FIX, no blocker or major), sourced at the end of dev/rel076-release-fix.R
# (local = TRUE there too, so rd(), wr(), edit_file() and restore() are the
# resolver's). The coordinator asked for findings 1 to 5 and 7 as release
# fixes R4 to R9; finding 6 (the ledger reasons of brmsfit-methods:923,
# :927, :937) stays as filed. R4 (the NEWS bullet) is in
# dev/rel076-news-split.R, which the NEWS merge applies after this file;
# R5's record and backlog text are in dev/round-20261027.md and
# dev/rel076-backlog-section.md; R9 is dev/rel076-commit.sh. The edits
# here change code, docs and tests. Each edit checks its old text is
# there once; a file no lane changed is restored from its HEAD blob first,
# so a merge rerun applies each edit once.
#
#   (sourced; not run on its own)

# ---- R5. dev/rtmb-pitfalls.md: the tmbstan pointer trap ----------------
# finding 2: dev/rel076-overlap-trip.R compared two tmbstan fits'
# log_prob() after both had run; tmbstan's sampling() keeps the objective
# in process-global pointers, so both calls evaluated the second run's
# objective
restore("dev/rtmb-pitfalls.md")
pf <- rd("dev/rtmb-pitfalls.md")
stopifnot(!any(grepl("^23[.] ", pf)), any(grepl("^22[.] ", pf)))
wr(c(pf, "",
"23. tmbstan holds the objective in process-global pointers: its",
"    `sampling()` method calls `.Call(\"set_pointers\", ...)` for the fit",
"    it runs. So `rstan::log_prob()` or `grad_log_prob()` of an EARLIER",
"    tmbstan fit in the same process evaluates the LATEST fit's",
"    objective, and a comparison of two fits made in one process",
"    compares one model with itself. At the 0.76.0 consolidation this",
"    made two runs' log densities look bitwise equal at 800 points;",
"    evaluated right after each run, 9 of 500 differ by one ulp",
"    (dev/reviews/2026-10-27-release.md, finding 2). Evaluate a",
"    stanfit's objective right after its own run, or in its own",
"    process."), "dev/rtmb-pitfalls.md")
cat("appended pitfall 23 to dev/rtmb-pitfalls.md\n")

# ---- R6. ac_bound_theta(): a bound at a computed end, within its ulps ---
# finding 3: 1 / max eigenvalue of a row-standardized W computes to
# 0.999999999999998 or to 1 by how LAPACK rounds, and ub = 1 was refused
# on exact inequality, printed as "(-1.72756, 1)". A bound within
# autocor_bound_chart()'s snap (64 ulps of the end) is that end, and a
# refusal prints the ends it compares with 15 digits.
edit_file("R/priors.R",
'  if (identical(what, "lb") && v == ends[1L]) return(-Inf)
  if (identical(what, "ub") && v == ends[2L]) return(Inf)',
'  # an end is computed (1 / an eigenvalue for sar()) and carries a few
  # ulps, so a bound typed at it, ub = 1 for a row-standardized W, is
  # that end within the snap autocor_bound_chart() applies
  at_end <- function(e) {
    abs(v - e) <= 64 * .Machine$double.eps * max(1, abs(e))
  }
  if (identical(what, "lb") && at_end(ends[1L])) return(-Inf)
  if (identical(what, "ub") && at_end(ends[2L])) return(Inf)')
edit_file("R/priors.R",
'               "the interval of ", ac[["label"]], ", (",
               format(tr$lo, digits = 6), ", ", format(tr$hi, digits = 6),
               ")", call. = FALSE)',
'               "the interval of ", ac[["label"]], ", (",
               format(tr$lo, digits = 15), ", ",
               format(tr$hi, digits = 15), ")", call. = FALSE)')
f <- "extensions/frmtmb.spatial/tests/testthat/test-sar.R"
x <- rd(f)
stopifnot(!any(grepl("a bound typed at a computed end", x, fixed = TRUE)))
wr(c(x, "", rd("dev/rel076-release-fix2-sar-test.R")), f)
cat("appended the computed-end test to", f, "\n")

# ---- R7. ?frm_sample: the sigma intercept under sigma ~ x --------------
# finding 4: R1's sentence in frmtmb.sample's help, the default included
edit_file("extensions/frmtmb.sample/R/sample.R",
"#' refused: a box on the innovation sd is not a box on the marginal sd,
#' which the sampler holds.",
"#' refused: a box on the innovation sd is not a box on the marginal sd,
#' which the sampler holds. Under `sigma ~ x` with such a term, a prior on
#' the intercept of `sigma` (`class = \"Intercept\", dpar = \"sigma\"`),
#' the default included, is still read at frmtmb's log marginal sd, where
#' brms reads it at the log innovation sd; the slopes are the same in
#' both, as [frmtmb::set_prior()] says.")

# ---- R8. log_lik(newdata =): name only an alternative that answers ------
# finding 5: the refusal said "Use posterior_epred(newdata =) or
# posterior_predict(newdata =)", and posterior_predict(newdata =) is
# refused for a model whose draws are structured (a correlated residual
# among them). The predicate posterior_predict() asks is one helper now,
# so the two methods cannot answer it apart. methods-draws.R is no
# lane's, so it is restored from its HEAD blob first.
md <- "extensions/frmtmb.sample/R/methods-draws.R"
restore(md)
edit_file(md,
'  structured <- local({
    ctx0 <- sim_context(fit, rspec, list(), aterms = av)
    if (arma_cond) ctx0[["autocor"]] <- NULL
    sim_is_structured(ctx0)
  })',
'  structured <- draws_pp_structured(fit, resp, av)')
edit_file(md,
"draws_refuse_newdata <- function(newdata, re_formula, re.form, what,",
"#' Whether `posterior_predict()` draws response `resp` as a walk over
#' the fitted rows (a hidden state sequence, a group-level latent
#' class, a correlated residual), and so refuses `newdata` and
#' `re_formula` for it. brms's `cov = FALSE` ARMA is not one: its rows
#' are drawn one by one around their one-step means.
#'
#' @noRd
draws_pp_structured <- function(fit, resp,
                                av = fit$frame[[\"aterm_values\"]][[resp]]) {
  ctx0 <- sim_context(fit, fit$spec$responses[[resp]], list(),
                      aterms = av)
  if (resp %in% arma_cond_resp(fit)) ctx0[[\"autocor\"]] <- NULL
  sim_is_structured(ctx0)
}

draws_refuse_newdata <- function(newdata, re_formula, re.form, what,")
edit_file("extensions/frmtmb.sample/R/loo.R",
'  draws_refuse_newdata(
    newdata, re_formula, arg_unset(), "log_lik()",',
'  # posterior_epred(newdata =) answers for every model, and
  # posterior_predict(newdata =) only where its draws are not structured
  alt <- paste0("Use posterior_epred(newdata =) or posterior_predict(",
                "newdata =) for a quantity on other rows")
  if (!is.null(newdata)) {
    fb <- draws_base_fit(object)
    if (any(vapply(names(fb$spec$responses), function(r) {
      draws_pp_structured(fb, r)
    }, NA))) {
      alt <- paste0("Use posterior_epred(newdata =) for a quantity on ",
                    "other rows; posterior_predict(newdata =) is ",
                    "refused for this model too, since its draws are ",
                    "structured")
    }
  }
  draws_refuse_newdata(
    newdata, re_formula, arg_unset(), "log_lik()",')
edit_file("extensions/frmtmb.sample/R/loo.R",
'    "re_formula would remove. Use posterior_epred(newdata =) or ",
    "posterior_predict(newdata =) for a quantity on other rows")',
'    "re_formula would remove. ", alt)')
f <- "extensions/frmtmb.sample/tests/testthat/test-loglik-autocor.R"
x <- rd(f)
stopifnot(!any(grepl("names only an alternative that answers", x,
                     fixed = TRUE)))
wr(c(x, "", rd("dev/rel076-release-fix2-loglik-test.R")), f)
cat("appended the log_lik(newdata =) test to", f, "\n")
