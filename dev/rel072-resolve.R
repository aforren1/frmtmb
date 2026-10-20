# Consolidation edits of the 0.72.0 release tree, applied after the
# textual merge (dev/rel072-merge.sh) and before the NEWS merge
# (dev/rel072-news-merge.R). Each edit checks that its old text is
# there exactly once, so a rerun on a fresh merge either applies all of
# them or stops on the first that no longer fits. A file no lane
# changed is restored from its HEAD blob (47cf659d, 0.71.0) before it
# is edited, so every edit here is applied once however often this
# runs. Adapted from dev/rel071-resolve.R.
#
#   Rscript dev/rel072-resolve.R
rel <- Sys.getenv("REL072_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-release")
# the checkout whose blobs are read (dev/rel072-merge.sh, REL072_GIT)
git_dir <- Sys.getenv("REL072_GIT", rel)
setwd(rel)
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)),
               stdout = TRUE)
  wr(sub("\r$", "", x), f)
}
edit_file <- function(f, old, new, n = 1L) {
  s <- paste(rd(f), collapse = "\n")
  hits <- gregexpr(old, s, fixed = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ",
                   substr(old, 1, 70))
  s <- gsub(old, new, s, fixed = TRUE)
  wr(strsplit(s, "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}
merged <- readLines("dev/rel072-log/merge-work/files")
restore <- function(f) {
  if (f %in% merged) stop(f, " is a merged file; edit it, do not restore it")
  from_head(f)
}

# ---- R/autoscale.R: illse2 against adhess, one line --------------------
# illse2 maps the joint precision back from the design's coordinates
# (sdr_unwhiten()); adhess hands sdreport() the Hessian in natural units
# (Hn, the AD Hessian when the check kept one). Both: the recipe of the
# adhess review (dev/reviews/2026-10-19-adhess.md, trial merges 1, 2).
edit_file("R/autoscale.R",
'<<<<<<< release
  sdr_rescue(fit, sdr_unwhiten(fit, sdr_name_joint(sdr, obj)),
             Hq / outer(u, u))
=======
  sdr_rescue(fit, sdr_name_joint(sdr, obj), Hn)
>>>>>>> adhess',
'  sdr_rescue(fit, sdr_unwhiten(fit, sdr_name_joint(sdr, obj)), Hn)')

# ---- R/confint.R: catmulti against drop1raw, drop1() -------------------
# Both paragraphs of the roxygen; in the code drop1raw's prior refusal
# comes first (it refuses the fit whatever its family), then catmulti's
# scope for a family without a mu, whose terms() replaces drop1raw's
edit_file("R/confint.R",
"<<<<<<< release
#' A family without a `mu` gives the main formula to every location",
"#' A family without a `mu` gives the main formula to every location")
edit_file("R/confint.R",
"#' drops it from every logit.
=======
#' Each refit reads the data the fit holds",
"#' drops it from every logit.
#'
#' Each refit reads the data the fit holds")
edit_file("R/confint.R",
"#' [anova.frmtmb_fit()].
>>>>>>> drop1raw",
"#' [anova.frmtmb_fit()].")
edit_file("R/confint.R",
"<<<<<<< release
  # the refits below edit the main formula, which a family without a mu",
"  if (fit_prior_density(object)) {
    prior_lrt_refusal(\"drop1()\", \"this fit has\")
  }
  # the refits below edit the main formula, which a family without a mu")
edit_file("R/confint.R",
"  tt <- terms(object, dpar = dp_main)
=======
  if (fit_prior_density(object)) {
    prior_lrt_refusal(\"drop1()\", \"this fit has\")
  }
  tt <- terms(object)
>>>>>>> drop1raw",
"  tt <- terms(object, dpar = dp_main)")

# ---- R/fit.R: the check_se memory sentence (adhess) --------------------
# adhess's NEWS and vignette say it; frmtmb_control()'s documentation is
# in illse2's file, so the lane left it to the consolidation
edit_file("R/fit.R",
"#'   parameters, where the check waited. `\"ignore\"` skips it. The
#'   warning names each parameter and the reason:",
"#'   parameters, where the check waited. `\"ignore\"` skips it. A model
#'   without random effects then builds the exact Hessian at the first
#'   standard-error use. That is faster than the finite differences it
#'   replaces but takes more memory: on a bernoulli model with 10,000
#'   rows the peak working set rose by 350 MB with 200 coefficients and
#'   by 870 MB with 500. The
#'   warning names each parameter and the reason:")

# ---- the reference-category twin after drop1raw (catmulti) -------------
# drop1raw's fit keeps the raw variables, so catmulti's reason for
# counting the twin at fit time (the model frame has no x for log(x))
# no longer holds; the count is still taken once, at the fit, so that
# logLik() never assembles a frame. The twin of a refit on the fit's
# own rows takes the fit's frozen bases and NA rows, as the refit does
edit_file("R/fit.R",
"  # a refcat = NA family's degrees of freedom are its reference twin's,
  # counted here while the raw data are in hand: the fit keeps only its
  # model frame, where log(x) has no x to re-evaluate (loglik_df())
  frame[[\"refcat_first_outer\"]] <- refcat_first_outer(
    spec, bform, data, na.action, data2, drop_unused_levels)",
"  # a refcat = NA family's degrees of freedom are its reference twin's,
  # counted once here, so that logLik() never assembles a frame
  # (loglik_df()); a refit on the fit's own rows gives the twin the
  # fit's bases and NA rows, as it gives the refit
  frame[[\"refcat_first_outer\"]] <- refcat_first_outer(
    spec, bform, data, na.action, data2, drop_unused_levels,
    predvar_map = refit[[\"predvar_map\"]],
    refit_na_action = refit[[\"na_action\"]])")
edit_file("R/methods-fit.R",
"#' frame from the same data while it has them (`refcat_first_outer()`)
#' and stores its parameter count; nothing here evaluates the formula
#' again, because the fit keeps only its model frame, where a term such
#' as `log(x)` has no raw `x` to read. Where the reference model is not",
"#' frame from the same data (`refcat_first_outer()`) and stores its
#' parameter count, so nothing here assembles a frame or evaluates the
#' formula again. Where the reference model is not")
edit_file("R/methods-fit.R",
"#' `loglik_df()`, from the data `frm()` was given. `NULL` for a model
#' without such a family; otherwise either the twin frame's template
#' lengths (`lens`, `nfix`) or the reason there is no count (`why`).
#'
#' @noRd
refcat_first_outer <- function(spec, bform, data, na.action, data2,
                               drop_unused_levels) {",
"#' `loglik_df()`, from the data `frm()` was given. `NULL` for a model
#' without such a family; otherwise either the twin frame's template
#' lengths (`lens`, `nfix`) or the reason there is no count (`why`).
#' `predvar_map` and `refit_na_action` are the refit channel's
#' (`refit_channel()`), passed to `assemble_frame()` as for the refit.
#'
#' @noRd
refcat_first_outer <- function(spec, bform, data, na.action, data2,
                               drop_unused_levels, predvar_map = NULL,
                               refit_na_action = NULL) {")
edit_file("R/methods-fit.R",
"    assemble_frame(parse_spec(bf2), data, na.action = na.action,
                   data2 = data2,
                   drop_unused_levels = drop_unused_levels),",
"    assemble_frame(parse_spec(bf2), data, na.action = na.action,
                   data2 = data2,
                   drop_unused_levels = drop_unused_levels,
                   predvar_map = predvar_map,
                   refit_na_action = refit_na_action),")

# ---- frmtmb.sample: the unreached-prior note (prilabel review R1) -------
# The shadowing prior can be the user's own (a coef row beside a class
# row), not only a default; the note names both
sp <- "extensions/frmtmb.sample/R/sample.R"
edit_file(sp,
'      spec_target(s), " applies to no parameter: a default written for ",
      "each of them is more specific and stays, as in brms. Name the ",
      "parameter\'s coef, as get_prior(route = \\"sample\\") lists it, to ",
      "replace that default"))',
'      spec_target(s), " applies to no parameter: a more specific prior ",
      "(a default or another specification) on each of them stays, as ",
      "in brms. To replace a default there, name the parameter\'s coef ",
      "as get_prior(route = \\"sample\\") lists it"))')
edit_file(sp,
"#' (dev/prilabel-rev-lscale.R). Such a specification applies nowhere,
#' which would pass in silence; the announcement names it and the coef
#' that reaches the parameter.",
"#' (dev/prilabel-rev-lscale.R). A user's own coef row shadows a
#' class-wide row of theirs the same way. Such a specification applies
#' nowhere, which would pass in silence; the announcement names it and
#' the coef that reaches the parameter.")
edit_file(sp,
"#' specification and stays, as in brms. The default-prior message
#' names a specification that applies to no parameter for this reason.",
"#' specification and stays, as in brms, and so does a more specific
#' specification of your own. The default-prior message names a
#' specification that applies to no parameter for either reason.")
edit_file("extensions/frmtmb.sample/tests/testthat/test-default-priors-brms.R",
'  # nor does a prior that reaches its parameters
  expect_false(grepl("applies to no parameter",
                     msgs(set_prior("normal(0, 1)", class = "sd",
                                    group = "g")), fixed = TRUE))
})',
'  # nor does a prior that reaches its parameters
  expect_false(grepl("applies to no parameter",
                     msgs(set_prior("normal(0, 1)", class = "sd",
                                    group = "g")), fixed = TRUE))
  # the release consolidation (lane prilabel\'s review, R1): a coef row
  # of the user\'s own shadows their class row too, and the note does
  # not blame a default for it
  m <- msgs(set_prior("normal(0, 1)", class = "sd") +
              set_prior("normal(0, 2)", class = "sd", group = "g",
                        coef = "Intercept"))
  expect_match(m, paste("class=sd applies to no parameter: a more",
                        "specific prior (a default or another",
                        "specification)"), fixed = TRUE)
})')

# ---- frmtmb's .onLoad(): TMB's finalizer (illse2 review, item 6) --------
restore("R/zzz.R")
edit_file("R/zzz.R",
".onLoad <- function(libname, pkgname) {
  frm_install_generics(pkgname)",
".onLoad <- function(libname, pkgname) {
  frm_install_generics(pkgname)
  # TMB's ADFun finalizer calls isNullPointer(), which TMB lazy-loads.
  # When a garbage collection runs the first finalizer of a session, it
  # forces that promise, and every finalizer pending with it fails with
  # \"promise already under evaluation\" (printed to the console) and
  # leaks its tape (dev/illse2-findings.md, dev/upstream-bugs.md TMB-4).
  # Forcing it here, before any tape exists, costs one lazy-load fetch
  tryCatch(invisible(get(\"isNullPointer\", envir = asNamespace(\"TMB\"),
                         inherits = FALSE)),
           error = function(e) NULL)")

# ---- pkgdown.yaml: the contributed compatibility rows (testiso R3) -----
restore(".github/workflows/pkgdown.yaml")
edit_file(".github/workflows/pkgdown.yaml",
"            if (st != 0L) stop(\"R CMD INSTALL failed for \", nm)
          }
        shell: Rscript {0}",
"            if (st != 0L) stop(\"R CMD INSTALL failed for \", nm)
          }
        shell: Rscript {0}
      # The rows the extensions contribute to the compatibility table
      # can be checked only with all seven loaded. No other job has them
      # all: coverage runs core alone, and R CMD check hides what the
      # DESCRIPTION does not name, so test-compat-extensions.R skips
      # there. This job has just installed them. A skip here would hide a
      # broken row, so it fails the step as a failure does
      - name: Check the extensions' compatibility rows
        env:
          NOT_CRAN: \"true\"
        run: |
          r <- as.data.frame(testthat::test_file(
            \"tests/testthat/test-compat-extensions.R\", package = \"frmtmb\",
            load_package = \"installed\", reporter = \"summary\"))
          print(r[, c(\"test\", \"nb\", \"failed\", \"skipped\", \"error\")])
          bad <- any(r$skipped) || any(r$error) || sum(r$failed) > 0 ||
            sum(r$nb) == 0
          if (bad) stop(\"test-compat-extensions.R did not pass with all \",
                        \"seven extensions installed\")
        shell: Rscript {0}")

# ---- versions and floors ------------------------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
restore("DESCRIPTION")
ver("DESCRIPTION", "0.71.0", "0.72.0")
restore("codemeta.json")
edit_file("codemeta.json", '"version": "0.71.0",', '"version": "0.72.0",')
for (p in c("frmtmb.sample", "frmtmb.eam", "frmtmb.latent",
            "frmtmb.coupling")) {
  restore(file.path("extensions", p, "DESCRIPTION"))
}
ds <- "extensions/frmtmb.sample/DESCRIPTION"
ver(ds, "0.19.0", "0.20.0")
edit_file(ds, "frmtmb (>= 0.71.0)", "frmtmb (>= 0.72.0)")
de <- "extensions/frmtmb.eam/DESCRIPTION"
ver(de, "0.12.0", "0.13.0")
edit_file(de, "frmtmb (>= 0.70.0)", "frmtmb (>= 0.72.0)")
# test-sample-ndt-default.R (lane sampleprior) needs frmtmb.sample's
# new ndt defaults
edit_file(de, "    frmtmb.sample,", "    frmtmb.sample (>= 0.20.0),")
dl <- "extensions/frmtmb.latent/DESCRIPTION"
ver(dl, "0.6.1", "0.6.2")
edit_file(dl, "frmtmb (>= 0.63.0)", "frmtmb (>= 0.72.0)")
dc <- "extensions/frmtmb.coupling/DESCRIPTION"
ver(dc, "0.6.0", "0.6.1")
# no floor: measured, its whole suite passes on frmtmb 0.71.0
# (dev/rel072-floor.sh, dev/round-20261020.md "Versions")


# ---- dev/test-backlog.md: closures in place, then this release's section
bl <- "dev/test-backlog.md"
restore(bl)
close_item <- function(old, note) edit_file(bl, old, paste0(old, note))
close_item("  sd. A `||` block could honor it. A priors lane.",
           " CLOSED at 0.72.0 (lane prilabel).")
close_item("  refused. Pre-existing for `|` terms (`dev/dblbar-preexist.R`).",
           " CLOSED at\n  0.72.0 (lane prilabel).")
close_item("  `?frmtmb-student-re` says to write one term to share the scale.",
           "\n  CLOSED at 0.72.0 (lane prilabel).")
close_item("  comes first, on rellib-r8 too (`dev/dblbar-rev2-seed2.R`).",
           " CLOSED at\n  0.72.0 (lane sampleprior).")
close_item("  `re_formula = ~ diag(x | g)` names no term. Pre-existing.",
           " CLOSED at\n  0.72.0 (lane prilabel), but for group-level `cs()`, which lane\n  prilabel filed again.")
close_item("  change is not cosmetic.",
           " The class \"cor\" refusal names the term as\n  written since 0.72.0 (lane prilabel); the tables still read it.")
close_item("  `variables(draws)` and `as_draws_df()`.",
           " CLOSED at 0.72.0 (lane\n  sampleprior: the draws carry brms's `sd_` and `cor_` columns).")
close_item("  multivariate frame); `dirichlet()` in `mvbf()` is untested.",
           " CLOSED at\n  0.72.0 (lane catmulti).")
close_item("  matrix for every multi-predictor family; brms gives `S x N x (K - 1)`.",
           "\n  CLOSED at 0.72.0 (lane catmulti).")
close_item("  object loses it.", " CLOSED at 0.72.0 (lane catmulti).")
close_item("- `multinomial()`'s `fitted()` refuses; brms returns `n x 4 x K`.",
           "\n  CLOSED at 0.72.0 (lane catmulti).")
close_item("  across families.", " CLOSED at 0.72.0 (lane sampleprior).")
close_item("  `marginaleffects::predictions()` on a matrix response fail.",
           paste0("\n  The `confint()` aliases and `drop1()` are CLOSED at 0.72.0",
                  " (lane\n  catmulti); `nobs()` is not a defect (brms counts",
                  " the same rows,\n  `dev/catmulti-findings.md`);",
                  " marginaleffects is filed again\n  (lane catmulti)."))
close_item("  fits, not 78. The coordinator filed it at the consolidation.",
           paste0(" CLOSED\n  at 0.72.0 (lane illse2's sd escape: the fit",
                  " ends 7.6e-9 short,\n  `dev/illse2-findings.md`,",
                  " section 2)."))
close_item("  read the whole matrix; at degrees 7 and 8 they stop with an error.",
           "\n  CLOSED at 0.72.0 (lane illse2).")
close_item("- REML, and nonlinear formulas, keep the \"illcond\" verdict.",
           " CLOSED at\n  0.72.0 (lane illse2).")
close_item("  seeds; on the others the route does not apply).",
           " The route reaches 400\n  outer parameters and bounds since 0.72.0 (lane illse2).")
close_item("  the factor.",
           " Not done at 0.72.0: lane sampleprior found the premise\n  wrong (`dev/sampleprior-findings.md`, item 4).")
close_item("  parameters, and 1.47 to 1.52 with `(1 | g)` (section 10, item 7).",
           "\n  At 0.72.0 the cost is 0.87 to 0.89 of 0.71.0's (lane illse2).")
close_item("  gen_extreme_value): the cheapest bin-2 rows to port next.",
           " CLOSED at\n  0.72.0 (lane bin2port).")
# a stray line of the 0.71.0 consolidation's section file
edit_file(bl, "\nNA\n  denominator (`dev/brmsport-bin2-summary.md`).\n", "\n")
edit_file(bl, "\n## Reference\n",
          paste0("\n", paste(rd("dev/rel072-backlog-section.md"),
                             collapse = "\n"), "\n\n## Reference\n"))

# ---- dev/upstream-bugs.md: brms-28, brms-29, TMB-4, brms-26's scope ----
ub <- "dev/upstream-bugs.md"
restore(ub)
edit_file(ub, "- Record: `dev/famsimplex-findings.md`, section 6, brms-B.",
          paste0("- Correction (lane catmulti, 2026-10-14;",
                 " `dev/catmulti-brms-b.R`, log\n  `brms-b.txt`): it",
                 " fires for any category with a formula of its own,\n",
                 "  not only the first predicted one. `mua ~ z` and",
                 " `mub ~ z` under\n  `refcat = \"c\"`, and `muc ~ z`",
                 " under the default reference, all stop\n  with the same",
                 " error. `dev/famsimplex-findings.md`, section 6, says\n",
                 "  \"the FIRST predicted category\".\n",
                 "- Record: `dev/famsimplex-findings.md`, section 6, brms-B;\n",
                 "  `dev/catmulti-findings.md`, section 13."))
edit_file(ub, "\n## RTMB\n",
          paste0("\n", paste(rd("dev/rel072-upstream.md"), collapse = "\n"),
                 "\n\n## RTMB\n"))
edit_file(ub, "\n## tmbstan\n",
          paste0("\n", paste(rd("dev/rel072-upstream-tmb.md"),
                             collapse = "\n"), "\n\n## tmbstan\n"))

# ---- the release harness reads rellib-r10 ------------------------------
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R")) {
  if (!f %in% merged) restore(f)
  edit_file(f, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/source/r/rellib-r10")
}
# testiso's new runner (copied untracked, so edited after each copy)
edit_file("dev/release/run-session.R", "C:/Users/adf44/source/r/rellib-r9",
          "C:/Users/adf44/source/r/rellib-r10")
for (f in c("dev/brmsport-bin2-all.sh", "dev/brmsport-bin2-record.sh",
            "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  wr(strsplit(gsub("rellib-r9", "rellib-r10", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}

# ---- bin 2: the verdict moves of this consolidation ---------------------
edit_file("dev/brmsport-bin2-verdicts-2p.R",
"W <- do.call(rbind, W)
stopifnot(!anyDuplicated(W$id)",
"W <- do.call(rbind, W)
# the verdict moves of the 0.72.0 consolidation
mv072 <- source(\"dev/rel072-bin2-moves.R\")$value
W <- mv072(W, \"values\")
V <- mv072(V, \"2p\")
stopifnot(!anyDuplicated(W$id)")

# ---- found at the consolidation: a lost correlation took the sds' -------
# intervals with it (lane adhess's sdr_nan_nonpos(), fix round 1, on
# test-review-v29.R's correlation run to 1; dev/rel072-v29-diag.R). Its
# NaN row and column spread through the delta method's product to
# every row of the block, the spread the comment above this code
# already names for a lost direction of the SE check
edit_file("R/confint.R",
'    pos <- th_pos[bk[["theta_idx"]]]
    se_g <- sqrt(pmax(diag(J %*% Vprop[pos, pos, drop = FALSE] %*% t(J)),
                      0))',
'    pos <- th_pos[bk[["theta_idx"]]]
    # A parameter whose variance is not shown (sdr_nan_nonpos(): a NaN
    # row and column, reason "nocurv") spreads its NaN through the
    # product as a lost direction does. A row that does not move along
    # it keeps its interval; a row that does gets none. The jacobian is
    # a central difference, so a row moved by a few ulps of its value
    # (rounding: the log sd of a us block against the parameter of its
    # correlation) has not moved
    Vb <- Vprop[pos, pos, drop = FALSE]
    nb <- which(!is.finite(diag(Vb)))
    kb <- setdiff(seq_along(pos), nb)
    Jk <- J[, kb, drop = FALSE]
    se_g <- sqrt(pmax(diag(Jk %*% Vb[kb, kb, drop = FALSE] %*% t(Jk)), 0))
    if (length(nb)) {
      hb <- 2e-5 * pmax(abs(t0[nb]), 1)
      tol <- 16 * .Machine$double.eps * pmax(abs(g0), 1)
      moves <- abs(sweep(J[, nb, drop = FALSE], 2, hb, `*`)) > tol
      se_g[rowSums(moves) > 0] <- NA_real_
    }')

# ---- dev/warnleak-scan.R takes FRMTMB_RELLIB, as run-tests.R does -------
# so that an interim build (wt-release-lib) is scanned before the
# release library exists
edit_file("dev/warnleak-scan.R",
'.libPaths(c("C:/Users/adf44/source/r/rellib-r10",',
'.libPaths(c(if (nzchar(Sys.getenv("FRMTMB_RELLIB"))) {
  Sys.getenv("FRMTMB_RELLIB")
}, "C:/Users/adf44/source/r/rellib-r10",')

# ---- found at the consolidation: sampleprior's draw-name test against ---
# prilabel's variables() order. Lane prilabel puts every block's cor_
# after every block's sd_, brms's variables() order
# (dev/prilabel-findings.md, R/confint.R hyp_env_vals(); brms-vars.txt);
# lane sampleprior's test wrote 0.71.0's interleaved order as brms's
edit_file("extensions/frmtmb.sample/tests/testthat/test-draws-sd-names.R",
'  # brms 2.23.0\'s names for this model
  expect_identical(want, c("sd_g__Intercept", "sd_g__x",
                           "cor_g__Intercept__x", "sd_h__Intercept"))',
'  # brms 2.23.0\'s names for this model, in its order: every sd_, then
  # every cor_ (lane prilabel, dev/prilabel-findings.md)
  expect_identical(want, c("sd_g__Intercept", "sd_g__x",
                           "sd_h__Intercept", "cor_g__Intercept__x"))')

# ---- dev/round-handoff.md: the 0.72.0 section (lane testiso changed the
# file, so it is merged and edited here)
hf <- "dev/round-handoff.md"
edit_file(hf,
"0.70.0 release and 2026-10-14 at the 0.71.0 release. Read this, then",
"0.70.0 release, 2026-10-14 at the 0.71.0 release and 2026-10-20 at
the 0.72.0 release. Read this, then")
edit_file(hf, "\n## The 0.71.0 round, 2026-10-13 to 2026-10-14\n",
          paste0("\n", paste(rd("dev/rel072-handoff-section.md"),
                             collapse = "\n"),
                 "\n\n## The 0.71.0 round, 2026-10-13 to 2026-10-14\n"))

# ---- the release review (dev/reviews/2026-10-20-release.md) ------------
# finding 5: lane catmulti assigned `multi` a second time, with the same
# value, in default_priors_for(); lane sampleprior's line at the top of
# the function stays
edit_file("extensions/frmtmb.sample/R/sample.R",
"  # own name (core's extra_map) and names its response, as brms does
  multi <- length(fit$spec$responses) > 1L
  for (r in names(fit$spec$responses)) {",
"  # own name (core's extra_map) and names its response, as brms does
  for (r in names(fit$spec$responses)) {")
# finding 4: US English
edit_file("R/predict.R", "#' `dec()` factor relevelled, or given",
          "#' `dec()` factor releveled, or given")
edit_file("tests/testthat/test-newdata-aterms.R",
          "  # 0.71.0: the relevelled rows swapped",
          "  # 0.71.0: the releveled rows swapped")
# finding 1: the consolidation's merged test is kept in dev/ and copied
# into place, so a fresh merge writes it
wr(rd("dev/rel072-test-refit-refcat.R"),
   "tests/testthat/test-refit-refcat.R")
cat("wrote tests/testthat/test-refit-refcat.R\n")
cat("RESOLVE DONE\n")
