# The hand resolutions of lane nanse's merge onto the 0.68.0 release
# tree, steps 3 to 5 and 8 to 14 of the nanse review's "Final merge
# recipe" (dev/reviews/2026-10-06-nanse.md in lane nanse's worktree).
# dev/rel068-nanse-merge.sh runs it after git merge-file has left the
# conflict blocks in the release files.
#
#   Rscript dev/rel068-nanse-resolve.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
rd <- function(f) readLines(file.path(root, f), warn = FALSE)
wr <- function(f, x) writeLines(x, file.path(root, f))

# conflict blocks as list(start, mid, end) line numbers
blocks <- function(x) {
  s <- grep("^<<<<<<< ", x)
  m <- grep("^=======$", x)
  e <- grep("^>>>>>>> ", x)
  stopifnot(length(s) == length(e))
  lapply(seq_along(s), function(i) {
    c(s[i], m[m > s[i] & m < e[i]][1], e[i])
  })
}
# replace block k of file f with `new` (a character vector); `ours`
# and `theirs` are passed to `new` when it is a function
resolve <- function(f, fns) {
  x <- rd(f)
  b <- blocks(x)
  stopifnot(length(b) == length(fns))
  for (k in rev(seq_along(b))) {
    bk <- b[[k]]
    ours <- if (bk[2] > bk[1] + 1L) x[(bk[1] + 1L):(bk[2] - 1L)] else
      character(0)
    theirs <- if (bk[3] > bk[2] + 1L) x[(bk[2] + 1L):(bk[3] - 1L)] else
      character(0)
    new <- fns[[k]](ours, theirs)
    x <- c(x[seq_len(bk[1] - 1L)], new, x[(bk[3] + 1L):length(x)])
  }
  wr(f, x)
  cat("resolved", length(b), "block(s) in", f, "\n")
}
both <- function(ours, theirs) c(ours, theirs)
edit <- function(f, old, new, all = FALSE) {
  x <- paste(rd(f), collapse = "\n")
  if (!grepl(old, x, fixed = TRUE)) stop("not found in ", f, ": ",
                                         substr(old, 1, 70))
  x <- if (all) gsub(old, new, x, fixed = TRUE) else
    sub(old, new, x, fixed = TRUE)
  wr(f, strsplit(x, "\n", fixed = TRUE)[[1]])
}

# 3. R/fit.R: the release's autoscaled, then the lane's start
resolve("R/fit.R", list(both))
# 4. R/interop.R: the release's extra_cov loop, then the lane's
#    jc_nonest block
resolve("R/interop.R", list(both))
# 5. R/predict.R: the release's v with ev_mu, then the lane's lost_k
#    lines, without the lane's own SE[, k] line
resolve("R/predict.R", list(function(ours, theirs) {
  c(ours, theirs[!grepl("^    SE\\[, k\\] <- ", theirs)])
}))
# 8. frmtmb.spline curve-cov.R, three hunks
resolve("extensions/frmtmb.spline/R/curve-cov.R", list(
  function(ours, theirs) {
    c(theirs[grepl("^  # ", theirs)],
      "  lost <- lb$se_nonest %||% rep(FALSE, nrow(C))",
      "  se <- unname(sqrt(pmax(diag(Sigma) + lb$extra_var, 0)))",
      "  se[lost] <- NaN",
      "  E <- lb$extra_cov",
      "  Sigma <- sp_add_extra(Sigma, E)",
      "  list(lb = lb, C = C, E = E, Sigma = Sigma, se = se, lost = lost,",
      "       span = lbc$span)")
  },
  function(ours, theirs) character(0),
  function(ours, theirs) {
    cm <- theirs[grepl("^  # ", theirs)]
    c(cm,
      "  lost <- st$lost[i1] | st$lost[i2]",
      "  se <- unname(sqrt(pmax(diag(Sigma), 0)))",
      "  se[lost] <- NaN",
      "  out <- list(eta = st$lb$eta[i1] - st$lb$eta[i2], C = C, V = st$lb$V,",
      "              Sigma = Sigma, E = E, se = se, lost = lost, rel = NA_real_,",
      "              n_predict = 0L, newdata = newdata, contrast = contrast,",
      "              dpar = dpar, resp = resp, re_formula = re_formula,",
      "              allow_new_levels = allow_new_levels, fit = fit, span = span,",
      "              extra_var = as.numeric(E[cbind(i1, i1)]))",
      "  if (!nl) {",
      "    r <- c(",
      "      sp_cov_check(fit, newdata, st$se[i1], dpar, resp, re_formula,",
      "                   allow_new_levels, tol, \"`newdata`\", st$lost[i1]),",
      "      sp_cov_check(fit, contrast, st$se[i2], dpar, resp, re_formula,",
      "                   allow_new_levels, tol, \"`contrast`\", st$lost[i2]))",
      "    out$rel <- if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)")
  }))
edit("extensions/frmtmb.spline/R/curve-cov.R",
     "                E = a$E, se = a$se, rel = NA_real_,",
     "                E = a$E, se = a$se, lost = a$lost, rel = NA_real_,")
# 9. curve-feature.R
resolve("extensions/frmtmb.spline/R/curve-feature.R", list(
  function(ours, theirs) {
    c("    .value = f0, .value_se = ifelse(lost, NaN, value_se),",
      "    stringsAsFactors = FALSE)")
  }))
# 10. test-curve.R: both tests, the release's first, closed
resolve("extensions/frmtmb.spline/tests/testthat/test-curve.R", list(
  function(ours, theirs) c(ours, "})", "", theirs)))

# 11. lane fixes' nl_flat_message(), the RB3 variant: flatness by
#     nl_flat_tol alone, the noise of the differenced Hessian in the
#     naming threshold only
edit("R/fit.R",
     "      H[, j] <- (obj$gr(up)[pos] - obj$gr(dn)[pos]) / (2 * h)\n    }\n    (H + t(H)) / 2\n  }, error = function(e) NULL)\n  if (is.null(H) || !all(is.finite(H))) return(NULL)",
     "      H[, j] <- (obj$gr(up)[pos] - obj$gr(dn)[pos]) / (2 * h)\n    }\n    H\n  }, error = function(e) NULL)\n  if (is.null(H) || !all(is.finite(H))) return(NULL)\n  # the asymmetry of the differenced matrix is its own noise (lane\n  # nanse's RB3): it sets how small a loading names a coefficient, and\n  # never whether a direction is flat\n  Hr <- H\n  H <- (Hr + t(Hr)) / 2")
edit("R/fit.R",
     "  S <- H[ok, ok, drop = FALSE] / outer(dg[ok], dg[ok])\n  ev <- eigen(S, symmetric = TRUE)\n  ratio <- abs(ev$values) / max(abs(ev$values))\n  flat <- ratio < nl_flat_tol\n  if (!any(flat)) return(NULL)\n  v <- ev$vectors[, flat, drop = FALSE]\n  load <- labels[ok][apply(abs(v) > 0.1, 1L, any)]",
     "  S <- H[ok, ok, drop = FALSE] / outer(dg[ok], dg[ok])\n  E <- abs(Hr - t(Hr))[ok, ok, drop = FALSE] / 2 / outer(dg[ok], dg[ok])\n  ev <- eigen(S, symmetric = TRUE)\n  big <- max(abs(ev$values))\n  flat <- abs(ev$values) < nl_flat_tol * big\n  if (!any(flat) || all(flat)) return(NULL)\n  # a coefficient is named by its projection on the flat subspace, so a\n  # direction shared by many coefficients names all of them\n  gap <- min(abs(ev$values[!flat]))\n  tau <- max(sqrt(nl_flat_tol), 10 * sqrt(sum(E^2)) / gap)\n  proj <- sqrt(rowSums(ev$vectors[, flat, drop = FALSE]^2))\n  load <- labels[ok][proj > tau]
  # no projection clears the noise of the differenced Hessian, so this
  # check cannot say which coefficients are flat; it stays silent and
  # the standard-error check names the parameters (the 0.68.0 merge:
  # the ridge a + b + log(c0) printed \"The coefficients  of the
  # nonlinear parameter '' ...\")
  if (!length(load)) return(NULL)")
# 12. the flat warning explains the lost standard errors, so lane
#     nanse's check does not repeat it
edit("R/fit.R",
     "    if (!is.null(flat_msg)) frm_warning(flat_msg, call. = FALSE)",
     "    if (!is.null(flat_msg)) {\n      frm_warning(flat_msg, call. = FALSE)\n      # explains every lost parameter but a bound-held one (se_report())\n      fit$cache$se_explained <- \"nl_flat\"\n    }")

# 13. tests: the fits where lane nanse's check now speaks (the edits of
#     the nanse review's trial merge, dev/nanse-rev3-merge-check/tree)
edit("tests/testthat/test-gp-by.R",
     "  fs <- frm(bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)), data = d)",
     paste0("  # two by-level GP sds run to 0 in sigma, so they have no standard\n",
            "  # error, and frm() says so; the test is about the joint density\n",
            "  lost <- \"Standard errors are not available\"\n",
            "  fs <- allow_warnings(frm(bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)),\n",
            "                           data = d), lost)"))
edit("tests/testthat/test-gp-by.R",
     paste0("  fb <- frm(mvbf(bf(y ~ gp(x, by = f, k = 8)),\n",
            "                 bf(y2 ~ gp(x, by = f, k = 6))) + set_rescor(FALSE),\n",
            "            data = d, family = gaussian())"),
     paste0("  fb <- allow_warnings(frm(mvbf(bf(y ~ gp(x, by = f, k = 8)),\n",
            "                                bf(y2 ~ gp(x, by = f, k = 6))) +\n",
            "                             set_rescor(FALSE),\n",
            "                           data = d, family = gaussian()), lost)"))
edit("tests/testthat/test-ordinal-mixture.R",
     "                          c(\"NaN\", \"Hessian\", \"gradient\", \"converge\"))",
     paste0("                          c(\"NaN\", \"Hessian\", \"gradient\", \"converge\",\n",
            "                            \"Standard errors are not available\"))"))
edit("extensions/frmtmb.spline/tests/testthat/test-difference.R",
     paste0("  fit <- frmtmb::frm(frmtmb::bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),\n",
            "                     family = stats::gaussian(), data = d)"),
     paste0("  # the smoothing sds run to 0, so they have no standard error, and\n",
            "  # frm() says so; the difference is what is tested\n",
            "  fit <- withCallingHandlers(\n",
            "    frmtmb::frm(frmtmb::bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),\n",
            "                family = stats::gaussian(), data = d),\n",
            "    warning = function(w) {\n",
            "      if (grepl(\"Standard errors are not available\", conditionMessage(w),\n",
            "                fixed = TRUE)) invokeRestart(\"muffleWarning\")\n",
            "    })"))
edit("tests/testthat/test-se-check.R",
     "  hit <- grep(se_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)\n  expect_length(hit, 1L)\n  expect_match(hit, \"61 of 62 parameters\", fixed = TRUE)",
     paste0("  # on the merged tree lane fixes' flat-direction warning speaks first and\n",
            "  # explains the lost SEs; it must name all 60 a_f* coefficients\n",
            "  hit <- grep(\"are not identified: at the optimum\", r$warnings,\n",
            "              fixed = TRUE, value = TRUE)\n",
            "  expect_length(hit, 1L)\n",
            "  expect_length(grep(se_lost_phrase, r$warnings, fixed = TRUE), 0L)\n",
            "  expect_length(regmatches(hit, gregexpr(\"a_f[0-9]+\", hit))[[1]], 60L)"))

# 14. frmtmb.learn's tests require 0.68.0's check (test-engine.R), so
#     its floor moves and it takes a patch bump; frmtmb_control() has a
#     fourteenth field, check_se
edit("extensions/frmtmb.learn/DESCRIPTION", "    frmtmb (>= 0.63.0)",
     "    frmtmb (>= 0.68.0)")
edit("extensions/frmtmb.learn/DESCRIPTION", "Version: 0.7.1",
     "Version: 0.7.2")
x <- rd("extensions/frmtmb.learn/NEWS.md")
if (!identical(x[1], "# frmtmb.learn 0.7.2")) {
  wr("extensions/frmtmb.learn/NEWS.md", c(
    "# frmtmb.learn 0.7.2",
    "",
    "No change to the package's code. Needs frmtmb 0.68.0 for its tests.",
    "",
    "* `test-engine.R`, `test-reference.R` and `test-stan-identity.R` allow",
    "  frmtmb 0.68.0's warning that standard errors are not available,",
    "  and `test-engine.R` requires it on its fit with a constant",
    "  likelihood, where no parameter has a standard error (lane nanse,",
    "  `dev/nanse-findings.md`).",
    "",
    x))
}
edit("extensions/frmtmb.sample/R/sample.R",
     "fifteen names above share none of frmtmb_control()'s thirteen.",
     "fifteen names above share none of frmtmb_control()'s fourteen.")
cat("RESOLVE DONE\n")
