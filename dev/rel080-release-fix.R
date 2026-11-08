# The merge steps of the brief that are edits, after dev/rel080-merge.sh:
#
#   Rscript dev/rel080-release-fix.R
#
# Each edit replaces text that must occur exactly once, and stops if it
# does not, so that a rerun on an edited tree stops at its first edit.
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
rd <- function(p) paste(readLines(file.path(REL, p), warn = FALSE),
                        collapse = "\n")
wr <- function(x, p) {
  con <- file(file.path(REL, p), open = "wb")
  on.exit(close(con))
  writeLines(strsplit(x, "\n", fixed = TRUE)[[1L]], con, sep = "\r\n")
}
edit <- function(p, old, new, what) {
  x <- rd(p)
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(p, ": ", what, ": the old text occurs ", n, " times")
  x <- sub(old, new, x, fixed = TRUE)
  if (grepl("[.]R$", p)) invisible(parse(text = x, keep.source = FALSE))
  wr(x, p)
  cat("EDIT", p, "|", what, "\n")
}

# step 3: loglik x nlstart. The sign probe compares what the fits
# maximized, as nl_restart_round() and nl_fit_off_saddle() do; logLik()
# stops under REML with a prior on the fixed effects, and the probe must
# not.
edit("R/par-template.R",
"  ll0 <- as.numeric(stats::logLik(fit))
  ll1 <- suppressWarnings(as.numeric(stats::logLik(r)))",
"  ll0 <- as.numeric(logLik_penalized(fit))
  ll1 <- suppressWarnings(as.numeric(logLik_penalized(r)))",
"step 3, nl_sign_probe() reads logLik_penalized()")
edit("R/par-template.R",
"  ll <- nl_ll_pair(pr$ll0, pr$ll1)
  frm_message(
    \"The fit of the nonlinear model from its default start reached \",
    \"log-likelihood \", ll[1L], \". Of \", pr$n, \" starts with the signs of the \",
    \"body's coefficients varied, the best reached \", ll[2L],",
"  ll <- nl_ll_pair(pr$ll0, pr$ll1)
  # the values are what the fits maximized, which under a prior with a
  # density is not logLik()
  llw <- if (fit_penalized(pr$fit)) \"penalized log-likelihood\" else
    \"log-likelihood\"
  frm_message(
    \"The fit of the nonlinear model from its default start reached \",
    llw, \" \", ll[1L], \". Of \", pr$n, \" starts with the signs of the \",
    \"body's coefficients varied, the best reached \", llw, \" \", ll[2L],",
"step 3, the multistart message names the penalized value")
# the unsettled warning compares the same two values (the tape's is the
# objective, which holds the log prior)
edit("R/par-template.R",
"      \"pattern's run stopped: it ended at log-likelihood \", lt[2L],",
"      \"pattern's run stopped: it ended at \", llw, \" \", lt[2L],",
"step 3, the unsettled warning names the same value")

# step 6: bigx x nestquad. The nested rounds run on tapes of their own in
# natural coordinates, with a box around each anchor; the final tape
# takes the design's basis for escape_stationary() and mo_search().
edit("R/fit.R",
"    qf <- quad_nest_fit(nll, template, random, frame[[\"map\"]], lap_obj,
                        ctl_opt, bounds, par_units, frame, qnest, vb)
    obj <- qf$obj
    opt <- qf$opt",
"    qf <- quad_nest_fit(nll, template, random, frame[[\"map\"]], lap_obj,
                        ctl_opt, bounds, par_units, frame, qnest, vb)
    obj <- qf$obj
    opt <- qf$opt
    # each round's tape is a rule anchored at one point and is optimized
    # inside a box around it, in natural coordinates (a box is not a box
    # in the basis); the final tape takes the basis, so that the escapes
    # below step in it as the Laplace fit does
    if (!is.null(ob)) obj$env$.frmtmb_opt_basis <- ob",
"step 6, the basis on the final nested tape")

# step 8: mfxse x readers. At a new level the mean integrates the level's
# effect over its block's covariance, so it reads that block's theta.
edit("R/interop.R",
"  keep_re <- re_form_keeps(dots[[\"re_formula\"]])
  b_idx <- integer(0)",
"  keep_re <- re_form_keeps(dots[[\"re_formula\"]])
  # a row at a level the fit did not see reads the mean over that level's
  # effect (fitted_me_moments()), which reads its block's theta
  new_lv <- isTRUE(dots[[\"allow_new_levels\"]])
  b_idx <- integer(0)",
"step 8, interop_re_part() keys new levels")
edit("R/interop.R",
"    if (block_reads_theta(bk)) th_idx <- c(th_idx, bk[[\"theta_idx\"]])
  }",
"    if (block_reads_theta(bk) ||
          (new_lv && !cs %in% c(\"smooth\", \"gp\", \"hsgp\"))) {
      th_idx <- c(th_idx, bk[[\"theta_idx\"]])
    }
  }",
"step 8, the theta of every non-field block at a new level")

# step 13: the release runners read this release's library
for (p in c("dev/release/run-tests.R", "dev/release/run-session.R")) {
  edit(p, "LIB <- \"C:/Users/adf44/source/r/rellib-r17\"",
       "LIB <- \"C:/Users/adf44/source/r/rellib-r19\"", "step 13, rellib-r19")
}

# versions, floors (steps 4 and 5) and roxygen2's version (step 14)
edit("DESCRIPTION", "Version: 0.79.1", "Version: 0.80.0", "frmtmb 0.80.0")
ver <- list(
  frmtmb.sample = c("0.26.0", "0.27.0"),
  frmtmb.latent = c("0.7.2", "0.7.3"),
  frmtmb.ode = c("0.7.2", "0.7.3"),
  frmtmb.survival = c("0.2.1", "0.2.2"))
for (pk in names(ver)) {
  edit(file.path("extensions", pk, "DESCRIPTION"),
       paste0("Version: ", ver[[pk]][1L]), paste0("Version: ", ver[[pk]][2L]),
       paste(pk, ver[[pk]][2L]))
}
for (pk in c("frmtmb.sample", "frmtmb.latent")) {
  edit(file.path("extensions", pk, "DESCRIPTION"), "frmtmb (>= 0.79.0)",
       "frmtmb (>= 0.80.0)", paste(pk, "floor 0.80.0"))
}
for (pk in list.files(file.path(REL, "extensions"))) {
  edit(file.path("extensions", pk, "DESCRIPTION"),
       "Config/roxygen2/version: 8.1.0", "Config/roxygen2/version: 8.1.1",
       paste(pk, "roxygen2 8.1.1"))
}
cat("RELEASE FIX DONE\n")
