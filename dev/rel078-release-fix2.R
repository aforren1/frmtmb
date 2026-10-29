# The fixes of the 0.78.0 release review (dev/reviews/2026-10-29-release.md,
# FIX, no blocker or major), sourced at the end of
# dev/rel078-release-fix.R (local = TRUE: rd(), wr(), edit_file(),
# restore(), merged and copied are the resolver's), so a merge rerun
# applies them. Finding 1 (the equiv review's evidence scripts) is in
# dev/rel078-merge.sh's skip list and dev/rel078-commit.sh's pathspec;
# finding 4 (a write-free dry run) in dev/rel078-commit.sh; findings 5
# and 8 in the record and dev/rel078-backlog-section.md.
#
# F2 (finding 2): prior_summary() of draws printed release fix R1's one
#    merged class row in place of the two rows the user wrote. The merged
#    specification now carries the two it came from (`shown`), and
#    prior_summary() lists those, as frmtmb.sample 0.24.0 and brms do; the
#    sampler reads the merged one. Passed back, the listed rows merge
#    again to the same prior.
# F3 (finding 3): bayes_R2() and loo_R2() refuse a response named twice,
#    with the message of draws_resp_set() (one helper, so the message
#    template stays unique). They keep brms's spelling only, which
#    test-brms-pins.R and test-loo-predict.R pin, so they do not take
#    draws_resp_set() itself, which also takes frmtmb's spelling.
# F6 (finding 6): dev/release/run-par.sh puts Rtools last on PATH.
# F7 (finding 7): the stale comment above bin2_pass_note.

sa <- "extensions/frmtmb.sample/R/sample.R"
md <- "extensions/frmtmb.sample/R/methods-draws.R"
lo <- "extensions/frmtmb.sample/R/loo.R"

# ---- F2: the rows as written, beside the merged slot --------------------
edit_file(sa,
"  resp1 <- names(fit$spec$responses)[1L]
  specs <- unclass(pl)
  rewritten <- logical(length(specs))",
"  resp1 <- names(fit$spec$responses)[1L]
  specs <- unclass(pl)
  # the specifications as written, for prior_summary() of a merged pair
  orig <- specs
  rewritten <- logical(length(specs))")
edit_file(sa,
"    if (is.na(a$lb %||% NA)) a$lb <- b$lb %||% NA
    if (is.na(a$ub %||% NA)) a$ub <- b$ub %||% NA
    specs[at] <- list(a)",
"    if (is.na(a$lb %||% NA)) a$lb <- b$lb %||% NA
    if (is.na(a$ub %||% NA)) a$ub <- b$ub %||% NA
    # the sampler reads the merged slot; prior_summary() lists the rows
    # the user wrote, both of them, as brms's does
    a$shown <- prior_spec_as_written(orig[sort(c(rw, cl))])
    specs[at] <- list(a)")
edit_file(sa,
"#' The classes of brms's proper defaults whose every position the",
"#' The specifications a merged slot came from, as they were written
#' (`refuse_coef_bounds()` keeps them on the merged one as `shown`): a
#' specification that itself came from a merge gives its own.
#'
#' @noRd
prior_spec_as_written <- function(specs) {
  out <- list()
  for (s in specs) {
    out <- c(out, if (is.null(s$shown)) list(s) else s$shown)
  }
  out
}

#' A prior list with each merged slot listed as the rows it came from,
#' for `prior_summary()` of draws: the rows the user wrote, both
#' spellings of a nonlinear intercept included, as brms's
#' `prior_summary()` keeps them (dev/reviews/2026-10-29-release.md,
#' finding 2). Passed back as `prior =`, they merge again to the slot the
#' sampler read.
#'
#' @noRd
prior_rows_as_written <- function(pl) {
  specs <- unclass(pl)
  if (!any(vapply(specs, function(s) !is.null(s$shown), NA))) return(pl)
  out <- list()
  for (s in specs) {
    if (is.null(s$shown)) {
      out <- c(out, list(s))
    } else {
      out <- c(out, s$shown)
    }
  }
  a <- attributes(pl)
  a$names <- NULL
  attributes(out) <- a
  out
}

#' The classes of brms's proper defaults whose every position the")
edit_file(md,
"    cat(\"No priors were used (flat improper priors on the outer \",
        \"parameters).\\n\", sep = \"\")
    return(invisible(NULL))
  }
  pl
}",
"    cat(\"No priors were used (flat improper priors on the outer \",
        \"parameters).\\n\", sep = \"\")
    return(invisible(NULL))
  }
  prior_rows_as_written(pl)
}")

# ---- F3: a response named twice, in the R2 functions ----------------------
edit_file(lo,
"  s <- match(as.character(resp), stan)
  if (anyNA(s)) draws_refuse_resp(paste(stan, collapse = \", \"))
  s
}",
"  s <- match(as.character(resp), stan)
  if (anyNA(s)) draws_refuse_resp(paste(stan, collapse = \", \"))
  # as draws_resp_set() refuses it for every other method: bayes_R2()
  # answered two identical rows, and loo_R2() two different ones, from
  # two runs of its bootstrap
  if (anyDuplicated(s)) draws_refuse_resp_twice(stan[s[anyDuplicated(s)]])
  s
}

#' The refusal of a `resp` that names one response twice, which brms
#' refuses too (on a rescor model it looks for the correlation of a
#' response with itself); `name` is the name as the caller spells it.
#'
#' @noRd
draws_refuse_resp_twice <- function(name) {
  frm_stop(\"`resp` names the response '\", name, \"' twice. Name each \",
           \"response once\", call. = FALSE)
}")
edit_file(md,
"    frm_stop(\"`resp` names the response '\", hit[anyDuplicated(hit)],
             \"' twice. Name each response once\", call. = FALSE)",
"    draws_refuse_resp_twice(hit[anyDuplicated(hit)])")

# ---- the tests, each seen to fail before the fix --------------------------
tp <- "extensions/frmtmb.sample/tests/testthat/test-prior-nlpar-intercept.R"
wr(c(rd(tp), rd("dev/rel078-fix2-r1summary-test.R")), tp)
cat("appended dev/rel078-fix2-r1summary-test.R to", tp, "\n")
tr <- "extensions/frmtmb.sample/tests/testthat/test-rescor-me.R"
wr(c(rd(tr), rd("dev/rel078-fix2-resp-test.R")), tr)
cat("appended dev/rel078-fix2-resp-test.R to", tr, "\n")

# ---- F6: dev/release/run-par.sh, Rtools last on PATH ----------------------
rp <- "dev/release/run-par.sh"
restore(rp)
edit_file(rp,
          "export PATH=\"$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH\"",
          "export PATH=\"$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin\"")

# ---- F7: the comment above bin2_pass_note ---------------------------------
# a comment of the spec; the tier files carry the note's text, not the
# comment, so they do not change (dev/rel078-merge.sh compares them)
edit_file("dev/brmsport-bin2-adapters-spec.R",
"# A literal pass that a reader could take for more than it shows. The
# two rescor posterior_predict rows are frmtmb's (core) rows, as every
# bin-2P row is; frmtmb.sample's posterior_predict() of set_rescor(TRUE)
# draws is another path, with a defect filed at the 0.77.0 release.",
"# A literal pass that a reader could take for more than it shows. The
# two rescor posterior_predict rows are frmtmb's (core) rows, as every
# bin-2P row is; frmtmb.sample's posterior_predict() of set_rescor(TRUE)
# draws reads the same joint draw since 0.25.0, which closed the defect
# filed at the 0.77.0 release (the note says so).")
