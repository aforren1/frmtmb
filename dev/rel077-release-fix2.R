# Release fixes after the 0.77.0 release review
# (dev/reviews/2026-10-28-release.md, verdict FIX), sourced at the end of
# dev/rel077-release-fix.R (local = TRUE there too, so rd(), wr(),
# edit_file(), restore() and append_test() are the resolver's). Each
# edit checks that its old text is there once; a file no lane changed is
# restored from its HEAD blob first, so a merge rerun applies each edit
# once. NEWS edits of fact are in dev/rel077-news-split.R; the record,
# the backlog section and dev/rel077-commit.sh carry the rest.
#
#   (sourced; not run on its own)
pf <- "R/priors.R"

# ---- F1 (finding 1, MAJOR). The route mark through `[` and subset() ---
# `[.data.frame` drops the `route` attribute whenever the call gives a
# column index, and subset() always gives one, so subset(pf, ...) and
# pf[, cols] of a route = "fit" table were read as written and sampled
# (flat) in place of brms's defaults (dev/review-rel077-c-run/route.txt,
# seed 11). A `[` method keeps the class and the route whenever the
# result is still a prior table (it has a `prior` column); a selection
# without one is a plain data frame. The plain-data-frame mark of
# release fix R1 is widened to any row of frmtmb's printed "(flat)"
# (brms's own tables store ""), so dropping the `source` column or
# setting it to "user" after as.data.frame() no longer hides the table
# (the review's two smaller gaps).
edit_file(pf,
"#' @export
as.data.frame.frmtmb_prior_rows <- function(x, row.names = NULL,",
"#' @export
`[.frmtmb_prior_rows` <- function(x, ...) {
  # a selection of rows or columns of a prior table is still that
  # table, so its route goes with it (frmtmb.sample refuses a route =
  # \"fit\" table however it was cut down); a selection that dropped the
  # `prior` column is no prior table and becomes a plain data frame
  route <- attr(x, \"route\")
  out <- NextMethod()
  if (is.data.frame(out)) {
    if (\"prior\" %in% names(out)) {
      attr(out, \"route\") <- route
      class(out) <- c(\"frmtmb_prior_rows\", \"data.frame\")
    } else {
      attr(out, \"route\") <- NULL
      class(out) <- \"data.frame\"
    }
  }
  out
}

#' @export
as.data.frame.frmtmb_prior_rows <- function(x, row.names = NULL,")
edit_file(pf,
"  # a table of default priors that does not say which route it
  # describes is the ambiguity `route` exists to remove, so the label
  # travels with the object rather than with the call that made it. A
  # ROW subset keeps it, correctly: fewer rows of a fit-route table are
  # still a fit-route table. A COLUMN subset does not, because
  # `[.data.frame` drops the attribute there, and a table missing the
  # `prior` column has no route to claim anyway. The NULL branch prints
  # such a table, and any assembled by other means, without a label.",
"  # a table of default priors that does not say which route it
  # describes is the ambiguity `route` exists to remove, so the label
  # travels with the object rather than with the call that made it: a
  # row or column selection or subset() keeps it (`[.frmtmb_prior_rows`),
  # since fewer rows or columns of a fit-route table are still one. The
  # NULL branch prints a table assembled by other means without a
  # label.")
edit_file(pf,
"  # frmtmb's own table, coerced: its route is gone, and its (flat) rows
  # may be frm()'s defaults; brms's tables store \"\" there, not \"(flat)\"
  if (is.null(route) && \"source\" %in% names(x) &&
      any(x[[\"prior\"]] %in% \"(flat)\" & x[[\"source\"]] %in% \"default\")) {
    route <- \"table\"
  }",
"  # frmtmb's own table, coerced: its route is gone, and its (flat) rows
  # may be frm()'s defaults. brms's tables store \"\" there, not
  # \"(flat)\", so frmtmb's printed spelling marks the table whatever its
  # `source` column says, or if it has none
  if (is.null(route) && any(x[[\"prior\"]] %in% \"(flat)\")) {
    route <- \"table\"
  }")
edit_file(pf,
"#' through `+`, `c()` and a column edit. The mark \"table\" is a plain
#' data frame of a frmtmb table, whose route as.data.frame() dropped:",
"#' through `+`, `c()`, a column edit, and a row or column selection or
#' `subset()` (`[.frmtmb_prior_rows`). The mark \"table\" is a plain
#' data frame holding frmtmb's printed `(flat)`, whose route
#' as.data.frame() dropped:")

# ---- F3 (finding 3). c() of a specification and a prior table ----------
# 0.76.0 accepted c(set_prior(), default_prior()); after reslot it was
# refused as "a plain data frame". brms 2.23.0's c() of a brmsprior and a
# prior table combines the rows in either order, and brm() then refuses
# a duplicated slot (dev/rel077-fix2-brms-c.R, fix2-brms-c.txt). So c()
# reads a prior table, or a brms prior, as as.brmsprior() reads it, with
# its route, and a table first in c() reaches the same method.
edit_file(pf,
"c.frmtmb_priorlist <- function(...) {
  parts <- list(...)
  for (p in parts) {",
"c.frmtmb_priorlist <- function(...) {
  # a default_prior() or validate_prior() table, or a brms prior, is read
  # as as.brmsprior() reads it, its route kept: brms's c() of a brmsprior
  # and a prior table combines the rows
  parts <- lapply(list(...), function(p) {
    if (inherits(p, c(\"frmtmb_prior_rows\", \"brmsprior\"))) {
      as.brmsprior(p)
    } else p
  })
  for (p in parts) {")
edit_file(pf,
"#' Print a prior specification
#'",
"#' @export
c.frmtmb_prior_rows <- function(...) {
  # a table first in c(): the same combination, as brms's c() of a prior
  # table and a brmsprior gives
  do.call(c.frmtmb_priorlist, list(...))
}

#' Print a prior specification
#'")

# ---- F5 (finding 5, nit). A stale comment on the named list -----------
edit_file(pf,
"#' through no fault of the caller. The legacy named list and the
#' `\"flat\"` string pass through untouched.",
"#' through no fault of the caller. The `\"flat\"` string passes through
#' untouched.")

# ---- F2 (finding 2), F6 (finding 6). ?frmtmb-marginaleffects -----------
# vcov_cluster() refuses a residual correlation term, so HC0, HC1, stata
# and a cluster formula are refused under ar(cov = TRUE) too
# (dev/review-rel077-b-run/x3.txt); and the 98-column line of release
# fix S2 is rewrapped
mi <- "R/interop.R"
edit_file(mi,
"#' marginaleffects maps to HC3) weight each row by its leverage, the
#' diagonal of a hat matrix, which frmtmb does not compute. `\"HAC\"`, `\"NeweyWest\"` and `\"kernHAC\"`
#' need a time order of the rows. `\"satterthwaite\"` and
#' `\"kenward-roger\"` are denominator degrees of freedom for a linear
#' mixed model fitted by REML. `vcov_cluster()` refuses a fit whose
#' random effects cross the clusters, and so `\"HC0\"` and `\"HC1\"`
#' refuse such a fit too.",
"#' marginaleffects maps to HC3) weight each row by its leverage, the
#' diagonal of a hat matrix, which frmtmb does not compute. `\"HAC\"`,
#' `\"NeweyWest\"` and `\"kernHAC\"` need a time order of the rows.
#' `\"satterthwaite\"` and `\"kenward-roger\"` are denominator degrees of
#' freedom for a linear mixed model fitted by REML. `vcov_cluster()`
#' refuses a fit whose random effects cross the clusters, and so
#' `\"HC0\"` and `\"HC1\"` refuse such a fit too. It also refuses a
#' residual correlation term (`ar()`, `ma()` or `arma()`, either form,
#' and the other residual covariance terms), whose density is joint over
#' each group: under one, `\"HC0\"`, `\"HC1\"`, `\"stata\"` and a cluster
#' formula are refused, `vcov = TRUE` gives the model-based standard
#' errors, and [frm_bootstrap()] is the route to a robust one.")

# ---- the tests of F1 and F3 ---------------------------------------------
append_test("extensions/frmtmb.sample/tests/testthat/test-flat-prior.R",
            "dev/rel077-fix2-route-test.R", "subset() and a column selection")
append_test("tests/testthat/test-prior-one-slot.R",
            "dev/rel077-fix2-c-test.R", "c() reads a prior table")
restore("tests/testthat/test-get-prior-route.R")
# core's own test pinned that a column selection drops the route; it now
# keeps it while the `prior` column is there
edit_file("tests/testthat/test-get-prior-route.R",
"  # a COLUMN subset keeps the class but loses the attribute, because
  # `[.data.frame` drops it. That is wanted rather than a `[` method
  # that carries it through: a selection that has dropped the `prior`
  # column has no route left to describe. This is also the only
  # ordinary way to reach the print method's unlabeled branch.
  cols <- gp[, c(\"class\", \"coef\")]
  expect_s3_class(cols, \"frmtmb_prior_rows\")
  expect_null(attr(cols, \"route\"))
  expect_false(grepl(\"route =\", utils::capture.output(print(cols))[1],
                     fixed = TRUE))",
"  # a COLUMN selection that keeps the `prior` column is still the
  # route's table, so it keeps the label (`[.frmtmb_prior_rows`; through
  # 0.77.0's merge `[.data.frame` dropped it there, and frmtmb.sample
  # sampled (flat) in place of brms's defaults); one without it is a
  # plain data frame
  pc <- gp[, c(\"prior\", \"class\", \"coef\")]
  expect_s3_class(pc, \"frmtmb_prior_rows\")
  expect_identical(attr(pc, \"route\"), \"fit\")
  expect_identical(attr(subset(gp, class != \"sd\"), \"route\"), \"fit\")
  cols <- gp[, c(\"class\", \"coef\")]
  expect_identical(class(cols), \"data.frame\")
  expect_null(attr(cols, \"route\"))
  # the print method's unlabeled branch: a table with no route
  nl <- gp
  attr(nl, \"route\") <- NULL
  expect_false(grepl(\"route =\", utils::capture.output(print(nl))[1],
                     fixed = TRUE))")
# the hand-built arm of reslot's test, as release fix R1 left it, drops
# the `source` column of the fit-route table, which the widened mark now
# catches (the review's first small gap); a table built by hand in brms's
# spelling, "" for (flat), is still read as written
edit_file("extensions/frmtmb.sample/tests/testthat/test-flat-prior.R",
"  # a table with no route (built by hand, so with no `source` column)
  # is read as written: its (flat) rows are brms's (flat)
  hand <- as.brmsprior(as.data.frame(p)[, setdiff(names(p), \"source\")])",
"  # a table with no route, built by hand in brms's spelling (\"\" for
  # (flat), no `source` column), is read as written: its empty rows are
  # brms's (flat)
  hd <- as.data.frame(p)[, setdiff(names(p), \"source\")]
  hd$prior[hd$prior == \"(flat)\"] <- \"\"
  hand <- as.brmsprior(hd)")


# `[.frmtmb_prior_rows` forwards i, j and drop to `[.data.frame` through
# NextMethod(), which the swallow detector of test-arg-refusal.R does not
# read as touching `...`; it joins the file's exemption list with its
# reason, as .get_sigma.frmtmb_fit did (the full core tier, r077-ref2:
# 1 failure there, "[.frmtmb_prior_rows")
edit_file("tests/testthat/test-arg-refusal.R",
"  # with its own caller's dots
  \".get_sigma.frmtmb_fit\"
)",
"  # with its own caller's dots
  \".get_sigma.frmtmb_fit\",
  # `[` hands its i, j and drop to `[.data.frame` through NextMethod()
  \"[.frmtmb_prior_rows\"
)")
