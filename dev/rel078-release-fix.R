# The release fixes of the 0.78.0 consolidation, sourced by
# dev/rel078-resolve.R (local = TRUE: rd(), wr(), edit_file() and
# restore() are the resolver's). Each edit checks that its old text is
# there exactly once, so a rerun on a fresh merge applies every fix or
# stops on the first that no longer fits.
#
# R1 (meme's review, short check of fix round 2, S3): on a ~ 1, with both
#    set_prior(class = "Intercept", nlpar = "a") and set_prior(class =
#    "b", nlpar = "a"), frm_sample() kept the later row; frm(),
#    rellib-r15 and brms let the Intercept-spelled density win in either
#    order (brms: the coef = "Intercept" row beats its class row). The
#    pair is now merged by specificity, as frm() resolves it (measured:
#    dev/rel078-r1-frm.R, dev/rel078-log/r1-frm.txt): the Intercept
#    spelling's density, else the class row's; on each side, its bound,
#    else the class row's.
# R2 (the coordinator): the bin-2 pass notes of posterior_predict:69 and
#    :72 said frmtmb.sample draws the responses independently; since this
#    release it draws them jointly through core's predict_rescor_draw().
# The cross-lane test (meme x rescorpp, mvbf() with rescor and me()) is
# a new frmtmb.sample test file, test-rescor-me.R, from
# dev/rel078-fix-rescorme-test.R (dev/rel078-overlap.R is its probe).

# ---- R1: frmtmb.sample R/sample.R, refuse_coef_bounds() ---------------
sa <- "extensions/frmtmb.sample/R/sample.R"
stopifnot(sa %in% merged)
edit_file(sa,
"#' than an intercept, a bound on it is a bound on one coefficient, and
#' it is refused as brms's spelling of it is.
#'
#' @noRd
refuse_coef_bounds <- function(fit, pl) {",
"#' than an intercept, a bound on it is a bound on one coefficient, and
#' it is refused as brms's spelling of it is.
#'
#' Written beside brms's class row of the same parameter, the spelling
#' is the more specific of the two, as brms's `coef = \"Intercept\"` row
#' is: the pair is one slot, resolved as frm() resolves it in either
#' order. Its density stands, or the class row's where it has none, and
#' on each side its bound stands, or the class row's where it has none.
#'
#' @noRd
refuse_coef_bounds <- function(fit, pl) {")
edit_file(sa,
"  resp1 <- names(fit$spec$responses)[1L]
  specs <- unclass(pl)
  for (i in seq_along(specs)) {",
"  resp1 <- names(fit$spec$responses)[1L]
  specs <- unclass(pl)
  rewritten <- logical(length(specs))
  for (i in seq_along(specs)) {")
edit_file(sa,
"        s$coef <- \"\"
        s$written <- NULL
        specs[[i]] <- s
        next",
"        s$coef <- \"\"
        s$written <- NULL
        specs[[i]] <- s
        rewritten[i] <- TRUE
        next")
edit_file(sa,
"  # a class row written both ways names one slot; the later one stands,
  # as a later specification of a slot does when two are stacked
  key <- vapply(specs, spec_target, \"\")
  keep <- !duplicated(key, fromLast = TRUE)",
"  # a class row written both ways names one slot. The later one stood,
  # so the order of two specifications decided the prior, where frm()
  # and brms let the intercept spelling win in either order
  # (dev/reviews/2026-10-28-meme.md, S3); merged here by specificity
  key <- vapply(specs, spec_target, \"\")
  for (k in unique(key[duplicated(key)])) {
    at <- which(key == k)
    rw <- at[rewritten[at]]
    cl <- at[!rewritten[at]]
    if (length(rw) != 1L || length(cl) != 1L) next
    a <- specs[[rw]]
    b <- specs[[cl]]
    if (is.null(a$dist)) {
      a$dist <- b$dist
      a$prior <- b$prior
    }
    if (is.na(a$lb %||% NA)) a$lb <- b$lb %||% NA
    if (is.na(a$ub %||% NA)) a$ub <- b$ub %||% NA
    specs[at] <- list(a)
  }
  keep <- !duplicated(key, fromLast = TRUE)")

# the test blocks: both orders on both routes, a fit's own prior then the
# call's, the bound from either row; the reviewer's case
# (dev/review-meme-sc2.R, section 3)
tp <- "extensions/frmtmb.sample/tests/testthat/test-prior-nlpar-intercept.R"
stopifnot(tp %in% copied)
wr(c(rd(tp), "", rd("dev/rel078-fix-r1-test.R")), tp)
cat("appended dev/rel078-fix-r1-test.R to", tp, "\n")

# ---- the cross-lane test: mvbf() with rescor and me() ------------------
to <- "extensions/frmtmb.sample/tests/testthat/test-rescor-me.R"
stopifnot(!file.exists(to) || !(to %in% c(merged, copied)))
wr(rd("dev/rel078-fix-rescorme-test.R"), to)
cat("wrote", to, "from dev/rel078-fix-rescorme-test.R\n")

# ---- R2: the bin-2 pass notes of posterior_predict:69 and :72 ----------
sf <- "dev/brmsport-bin2-adapters-spec.R"
restore(sf)
edit_file(sf,
"  \"pass through core's predict() joint draw\",
  \"(brms_adapt_mv_posterior_predict()), the route of every bin-2P row;\",
  \"it says nothing about\",
  \"frmtmb.sample's posterior_predict() of set_rescor(TRUE) draws, which\",
  \"draws the responses independently (dev/test-backlog.md, filed at the\",
  \"0.77.0 release)\"), 2L),",
"  \"pass through core's predict() joint draw\",
  \"(brms_adapt_mv_posterior_predict()), the route of every bin-2P row.\",
  \"frmtmb.sample's posterior_predict() of set_rescor(TRUE) draws reads\",
  \"the same joint draw through core's predict_rescor_draw() since\",
  \"frmtmb.sample 0.25.0, which its own tests check\",
  \"(test-rescor-draws.R)\"), 2L),")

# ---- the release review's fixes (dev/reviews/2026-10-29-release.md) -----
source("dev/rel078-release-fix2.R", local = TRUE)
