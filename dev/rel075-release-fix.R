# Release fixes of the 0.75.0 consolidation, from the release review
# (dev/reviews/2026-10-24-release.md, FIX, no blocker), as the
# coordinator asked for them. Sourced by dev/rel075-resolve.R after the
# conflict resolutions (local = TRUE, so edit_file() and rd() are the
# resolve step's), so a merge rerun applies them again. Each edit
# matches its old text exactly once. The record's section "After the
# release review" has the evidence.

# ---- 1. (major) emmeans: no head check inside emmeans's try() ------------
# smallfix put require_fit_packages() at the head of
# recover_data.frmtmb_fit(). emmeans calls that function inside try()
# and replaces whatever it raises with "Perhaps a 'data' or 'params'
# argument is needed", so ref_grid() and emmeans() of a fit with a
# missing part lost the classed refusal 0.74.0 gave (review finding 1).
# emm_basis.frmtmb_fit(), which emmeans calls outside try(), keeps its
# head check and refuses by the package's name. The consolidation's
# record kept the head check on a wrong premise; the coordinator
# overruled it
edit_file("R/interop.R",
"                                    re_formula = NA, epred = FALSE) {
  require_fit_packages(object)
  # emmeans calls this inside try() and replaces whatever it raises with",
"                                    re_formula = NA, epred = FALSE) {
  # No require_fit_packages() here, though the other methods have one:
  # emmeans would discard it, as it discards every error raised here.
  # emm_basis.frmtmb_fit(), which emmeans calls outside try(), refuses a
  # fit with a missing part by the package's name.
  #
  # emmeans calls this inside try() and replaces whatever it raises with")
# the test: ref_grid() and emmeans() join the methods every block of
# test-extension-home-mixed.R calls
edit_file("tests/testthat/test-extension-home-mixed.R",
"  ngrps = function(f) ngrps(f), variables = function(f) variables(f))",
"  ngrps = function(f) ngrps(f), variables = function(f) variables(f),
  # emmeans runs recover_data() inside try() and replaced a refusal
  # raised there with \"Perhaps a 'data' or 'params' argument is needed\"
  # (the 0.75.0 release review, finding 1)
  ref_grid = function(f) emmeans::ref_grid(f),
  emmeans = function(f) emmeans::emmeans(f, ~ x))")

# ---- 2. frmtmb_covstruct_info() of a residual correlation term ------------
# sarfcor's homes table lists sar() among the extension specials, so the
# reader sent "sar" to cs_entry(), which with frmtmb.spatial loaded
# raised frmtmb_extension_missing ("The sar() structure is not
# registered in this session"; review finding 2). A residual term is
# recognized first and named as what it is
edit_file("R/covstruct-seam.R",
"#'   A name nothing registered is refused. A name that core lists as an
#'   extension's (`\"car\"`, `\"spde\"`) loads that extension, which
#'   registers it, as a fit that uses the structure would.",
"#'   A name nothing registered is refused. A name that core lists as an
#'   extension's (`\"car\"`, `\"spde\"`) loads that extension, which
#'   registers it, as a fit that uses the structure would. A residual
#'   correlation term (`\"fcor\"`, frmtmb.spatial's `\"sar\"`) is no
#'   covariance structure: it is refused with a pointer to
#'   [frmtmb_autocor_pars()].")
edit_file("R/covstruct-seam.R",
"  e <- frmtmb_covstruct_registry$reg[[name]]
  if (is.null(e) && !is.null(extension_special_home(name))) {",
"  # a residual correlation term is registered through another seam, and
  # an extension's (sar()) is in the same table of homes as a structure
  if (name %in% c(ac_entry_names(), extension_autocor_specials)) {
    frm_stop(\"frmtmb_covstruct_info(): '\", name, \"' is a residual \",
             \"correlation term, registered with frmtmb_register_autocor(), \",
             \"not a covariance structure. frmtmb_autocor_pars(fit) lists \",
             \"the parameters of a fit's residual correlation terms\",
             call. = FALSE)
  }
  e <- frmtmb_covstruct_registry$reg[[name]]
  if (is.null(e) && !is.null(extension_special_home(name))) {")
tw <- "tests/testthat/test-autocor-seam.R"
stopifnot(!any(grepl("residual term is not a structure", rd(tw),
                     fixed = TRUE)))
wr(c(rd(tw), "",
"test_that(\"frmtmb_covstruct_info() says a residual term is not a structure\", {",
"  # the 0.75.0 release review, finding 2: \"sar\" was sent to the",
"  # structures' registry and refused as a missing package",
"  local_frmtmb_registries()",
"  register_myfix()",
"  for (nm in c(\"fcor\", \"sar\", \"myfix\")) {",
"    e <- tryCatch(frmtmb_covstruct_info(nm), error = identity)",
"    expect_s3_class(e, \"frmtmb_error\")",
"    expect_false(inherits(e, \"frmtmb_extension_missing\"), label = nm)",
"    expect_match(conditionMessage(e), paste0(\"'\", nm, \"' is a residual \",",
"                                             \"correlation term\"),",
"                 fixed = TRUE)",
"    expect_match(conditionMessage(e), \"frmtmb_autocor_pars(fit)\",",
"                 fixed = TRUE)",
"  }",
"  # the guard absent: a structure and an unknown name read as before",
"  expect_false(frmtmb_covstruct_info(\"ar1\")$field)",
"  expect_error(frmtmb_covstruct_info(\"nosuchname\"),",
"               \"no covariance structure 'nosuchname'\", fixed = TRUE)",
"})"), tw)
cat("appended a block to", tw, "\n")

# ---- 3. one sar() parameter, one prior slot ----------------------------
# default_prior() listed lagsar and its internal parameter
# (class = "theta", coef = "thetaac_1"), validate_prior() and frm() took
# a prior on each, and the fit applied both (review finding 3). For a
# residual term's parameter with `bounds` (sar()'s lagsar or errorsar)
# the internal row leaves the table and a prior on it is refused by the
# class's name. ar()'s same double slot (pre-existing since the ar()
# classes) is filed in dev/test-backlog.md and not changed here
edit_file("R/autocor-seam.R",
"#' The class names `set_prior()` reaches on a frame's registry terms:",
"#' The internal names (`thetaac_<k>`) of a frame's bounded registry
#' parameters, named by their class (`lagsar`): such a parameter has
#' one prior slot, its class, so the internal spelling is not listed
#' in the prior table and a prior on it is refused.
#'
#' @noRd
autocor_seam_bounded_coefs <- function(frame) {
  out <- character(0)
  for (ac in frame[[\"autocor\"]] %||% list()) {
    for (pr in autocor_seam_bounded(ac)) {
      out[paste0(\"thetaac_\", ac[[\"theta_idx\"]][pr$idx])] <- pr$name
    }
  }
  out
}

#' The class names `set_prior()` reaches on a frame's registry terms:")
edit_file("R/priors.R",
"  if (length(frame[[\"par_template\"]][[\"theta\"]] %||% numeric(0))) add(\"theta\")
  for (cmp in theta_components) {
    v <- frame[[\"par_template\"]][[cmp]] %||% numeric(0)
    if (!length(v)) next
    for (nm in par_template_names(v, cmp)) add(\"theta\", coef = nm)
  }",
"  if (length(frame[[\"par_template\"]][[\"theta\"]] %||% numeric(0))) add(\"theta\")
  # a registry term's bounded parameter is listed once, under its class
  # (lagsar) above, as brms lists it
  own_class <- autocor_seam_bounded_coefs(frame)
  for (cmp in theta_components) {
    v <- frame[[\"par_template\"]][[cmp]] %||% numeric(0)
    if (!length(v)) next
    for (nm in par_template_names(v, cmp)) {
      if (identical(cmp, \"thetaac\") && nm %in% names(own_class)) next
      add(\"theta\", coef = nm)
    }
  }")
edit_file("R/priors.R",
"      tgt <- theta_coef_target(frame, s$coef)
      cmp <- tgt$comp
      nms <- par_template_names(frame[[\"par_template\"]][[cmp]] %||% numeric(0),
                                cmp)",
"      tgt <- theta_coef_target(frame, s$coef)
      cmp <- tgt$comp
      nms <- par_template_names(frame[[\"par_template\"]][[cmp]] %||% numeric(0),
                                cmp)
      # one slot per parameter: a bounded registry parameter's internal
      # name would be a second prior on the parameter its class names
      own_class <- autocor_seam_bounded_coefs(frame)
      hit <- if (identical(cmp, \"thetaac\")) {
        intersect(nms[tgt$idx], names(own_class))
      } else character(0)
      if (length(hit)) {
        frm_stop(\"class = \\\"theta\\\", coef = \\\"\", hit[1L], \"\\\" is the \",
                 \"internal parameter of \", own_class[[hit[1L]]], \", which \",
                 \"has a class of its own, as in brms. Write the prior with \",
                 \"class = \\\"\", own_class[[hit[1L]]], \"\\\" (lb and ub narrow \",
                 \"its interval); one parameter takes one prior\",
                 call. = FALSE)
      }")

# ---- 6. a residual term's class named in the refusal of another --------
# class = "lagsar" on an error SAR said the class "names no
# distributional parameter of this model" and listed neither the
# model's errorsar nor any registered class (review finding 6)
edit_file("R/priors.R",
"    return(paste0(spelling, \" names no distributional parameter of \",
                  \"this model. It has \", paste(have, collapse = \", \"),
                  \". frmtmb's own classes are \",
                  paste(frmtmb_prior_classes, collapse = \", \")))",
"    ac_cls <- autocor_seam_classes(frame)
    return(paste0(spelling, \" names no distributional parameter of \",
                  \"this model. It has \", paste(have, collapse = \", \"),
                  if (length(ac_cls)) {
                    paste0(\". Its residual correlation terms take \",
                           paste0(\"class = \\\"\", ac_cls, \"\\\"\",
                                  collapse = \", \"))
                  },
                  \". frmtmb's own classes are \",
                  paste(frmtmb_prior_classes, collapse = \", \")))")

# the tests of 3 and 6, in frmtmb.spatial's test-sar.R (sar() is its)
ts <- "extensions/frmtmb.spatial/tests/testthat/test-sar.R"
stopifnot(!any(grepl("one sar() parameter has one prior slot", rd(ts),
                     fixed = TRUE)))
wr(c(rd(ts), "",
"test_that(\"one sar() parameter has one prior slot, under its brms class\", {",
"  # the 0.75.0 release review, finding 3: default_prior() listed lagsar",
"  # and its internal parameter, and a prior on each was applied",
"  s <- sar_data(k = 5)",
"  dd <- s$d",
"  d2 <- list(W = s$W)",
"  dp <- default_prior(ylag ~ x + sar(W), data = dd, data2 = d2)",
"  expect_true(\"lagsar\" %in% dp$class)",
"  expect_false(any(dp$class == \"theta\"))",
"  for (type in c(\"lag\", \"error\")) {",
"    fo <- stats::as.formula(paste0(\"y\", if (type == \"lag\") \"lag\" else",
"                                   \"err\", \" ~ x + sar(W, type = '\", type,",
"                                   \"')\"))",
"    e <- tryCatch(validate_prior(set_prior(\"normal(0, 0.5)\",",
"                                           class = \"theta\",",
"                                           coef = \"thetaac_1\"),",
"                                 fo, data = dd, data2 = d2),",
"                  error = identity)",
"    expect_s3_class(e, \"frmtmb_error\")",
"    expect_match(conditionMessage(e), paste0(\"internal parameter of \",",
"                                             type, \"sar\"), fixed = TRUE)",
"    expect_match(conditionMessage(e), paste0(\"class = \\\"\", type, \"sar\\\"\"),",
"                 fixed = TRUE)",
"  }",
"  # the guard absent: the class itself takes a prior",
"  vp <- validate_prior(set_prior(\"normal(0.3, 0.2)\", class = \"lagsar\"),",
"                       ylag ~ x + sar(W), data = dd, data2 = d2)",
"  expect_true(\"normal(0.3, 0.2)\" %in% vp$prior[vp$class == \"lagsar\"])",
"  # and core's ar() keeps its internal row (its own double slot is filed",
"  # for a later round)",
"  da <- data.frame(t = 1:30, y = stats::rnorm(30))",
"  dpa <- default_prior(y ~ ar(time = t), data = da)",
"  expect_true(any(dpa$class == \"theta\" & dpa$coef == \"thetaac_1\"))",
"})",
"",
"test_that(\"class = 'lagsar' on an error SAR names the model's errorsar\", {",
"  # the 0.75.0 release review, finding 6",
"  s <- sar_data(k = 5)",
"  e <- tryCatch(validate_prior(set_prior(\"normal(0, 0.5)\", class = \"lagsar\"),",
"                               yerr ~ x + sar(W, type = \"error\"),",
"                               data = s$d, data2 = list(W = s$W)),",
"                error = identity)",
"  expect_s3_class(e, \"frmtmb_error\")",
"  expect_match(conditionMessage(e), paste0(\"Its residual correlation \",",
"                                           \"terms take class = \\\"errorsar\\\"\"),",
"               fixed = TRUE)",
"  # the guard absent: a model with no residual term says nothing of one",
"  e2 <- tryCatch(validate_prior(set_prior(\"normal(0, 0.5)\", class = \"lagsar\"),",
"                                yerr ~ x, data = s$d), error = identity)",
"  expect_s3_class(e2, \"frmtmb_error\")",
"  expect_no_match(conditionMessage(e2), \"residual correlation terms take\",",
"                  fixed = TRUE)",
"})",
"",
"test_that(\"frmtmb_covstruct_info() says sar() is a residual term\", {",
"  # the 0.75.0 release review, finding 2, with frmtmb.spatial loaded",
"  e <- tryCatch(frmtmb_covstruct_info(\"sar\"), error = identity)",
"  expect_s3_class(e, \"frmtmb_error\")",
"  expect_false(inherits(e, \"frmtmb_extension_missing\"))",
"  expect_match(conditionMessage(e), \"'sar' is a residual correlation term\",",
"               fixed = TRUE)",
"  # the guard absent: car(), a structure of this package, answers",
"  expect_true(frmtmb_covstruct_info(\"car\")$field)",
"})"), ts)
cat("appended three blocks to", ts, "\n")

# ---- 4. frmtmb.survival: bhaz() in the README and the vignette ----------
# the package's main new feature for a brms user was in NEWS, ?cox and
# the package help only (review finding 4)
sv <- "extensions/frmtmb.survival/README.md"
restore(sv)
edit_file(sv,
"  `cox_baseline()` gives the fitted weights.",
"  `cox_baseline()` gives the fitted weights.
- `bhaz()` is brms's addition term for the baseline basis of `cox()`:
  `time | cens(out) + bhaz(df = 6) ~ x` with `family = cox()` is the
  model of `cox(df = 6)`.")
vg <- "extensions/frmtmb.survival/vignettes/royston-parmar.Rmd"
restore(vg)
edit_file(vg,
"`royston_parmar()` is that family, parameterized exactly as
`flexsurv::flexsurvspline()` parameterizes it.",
"`royston_parmar()` is that family, parameterized exactly as
`flexsurv::flexsurvspline()` parameterizes it. The Cox model is
`cox()`; the last section shows brms's spelling of its baseline.")
wr(c(rd(vg), "",
"## brms's spelling of a Cox baseline",
"",
"`cox()` writes its baseline hazard as an M-spline with simplex weights,",
"as brms does. brms sets that basis with the addition term `bhaz()`,",
"and frmtmb.survival takes the same spelling. The two fits below are",
"one model:",
"",
"```{r bhaz}",
"f_bhaz <- frm(bf(recyrs | cens(censored) + bhaz(df = 6) ~ group),",
"              family = cox(), data = bc)",
"f_cox <- frm(bf(recyrs | cens(censored) ~ group),",
"             family = cox(df = 6), data = bc)",
"c(bhaz = logLik(f_bhaz), cox = logLik(f_cox))",
"```"), vg)
cat("appended a section to", vg, "\n")

# ---- 5. the comments on the saved-fit step name both fits ------------------
# the consolidation added a car() and sar() fit to the file the step runs
# (review finding 5)
edit_file(".github/workflows/pkgdown.yaml",
"      # The same rule for the saved-fit refusal: core's
      # test-extension-home-mixed.R reads a cox() and car() fit where
      # frmtmb.survival and frmtmb.spatial cannot be used, which needs
      # both installed to save it, so it skips under R CMD check and in
      # the coverage jobs.",
"      # The same rule for the saved-fit refusal: core's
      # test-extension-home-mixed.R reads a cox() and car() fit where
      # frmtmb.survival and frmtmb.spatial cannot be used, and a car()
      # and sar() fit where frmtmb.spatial cannot be used. Saving them
      # needs both installed, so those blocks skip under R CMD check and
      # in the coverage jobs.")
edit_file(".github/scripts/saved-fit-check.R",
"# The file's last blocks save a cox() and car() fit with frmtmb.survival
# and frmtmb.spatial loaded, then read it in child processes where they
# are absent.",
"# Two of the file's blocks save a fit in a child process, a cox() and
# car() fit with frmtmb.survival and frmtmb.spatial loaded and a car()
# and sar() fit with frmtmb.spatial loaded, then read it in child
# processes where those packages are absent.")
cat("RELEASE FIX DONE\n")
