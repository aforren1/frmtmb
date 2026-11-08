# Semantic conflicts found on the merged build, after dev/rel080-release-fix.R:
#
#   Rscript dev/rel080-release-fix2.R
#
# Same rule: each old text occurs exactly once.
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

# fitedge x nestquad: check_convergence()'s stop warnings now carry their
# own classes (frmtmb_conv_code, frmtmb_conv_gradient,
# frmtmb_conv_gradient_nonfinite) before frmtmb_warning, so the nested
# test's first class is no longer frmtmb_warning. It reads the class it
# means: the generic convergence warning, by inheritance.
edit("tests/testthat/test-quadrature-nested.R",
"    warning = function(w) {
      seen <<- c(seen, class(w)[1L])
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart(\"muffleWarning\")
    })
  # either the rule settled, or a warning of its own class says it did
  # not; a stop at the round cap with nothing said is the defect
  expect_true(isTRUE(fq$quad_nested$settled) ||
                \"frmtmb_quad_unsettled\" %in% seen)
  # every other warning here is the generic convergence one
  expect_true(all(seen %in% c(\"frmtmb_quad_unsettled\", \"frmtmb_warning\")))",
"    warning = function(w) {
      seen <<- c(seen, if (inherits(w, \"frmtmb_quad_unsettled\")) {
        \"frmtmb_quad_unsettled\"
      } else if (inherits(w, \"frmtmb_warning\")) \"frmtmb_warning\" else {
        class(w)[1L]
      })
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart(\"muffleWarning\")
    })
  # either the rule settled, or a warning of its own class says it did
  # not; a stop at the round cap with nothing said is the defect
  expect_true(isTRUE(fq$quad_nested$settled) ||
                \"frmtmb_quad_unsettled\" %in% seen)
  # every other warning here is a convergence warning of the package
  # (frmtmb_conv_code and frmtmb_conv_gradient are frmtmb_warning too)
  expect_true(all(seen %in% c(\"frmtmb_quad_unsettled\", \"frmtmb_warning\")))",
"fitedge x nestquad, the nested test reads frmtmb_warning by inheritance")
cat("RELEASE FIX 2 DONE\n")

# the user's decision of 2026-11-06 on TMB's internals (via the
# coordinator): keep both reads, and file the upstream request
edit("dev/test-backlog.md",
"### From the reviews, filed at the release\n\n",
paste0("### From the reviews, filed at the release\n\n",
"- **UPSTREAM, for the user. Ask TMB upstream to export `solveSubset()`\n",
"  and a tape-size accessor (`info()`); user will raise it later.** 0.80.0\n",
"  keeps both reads: `solveSubset()` through `get0()` with a correct\n",
"  fallback (latent residuals' selected inverse), `info()` through\n",
"  `getFromNamespace()` (nested quadrature's memory guard). Nothing is\n",
"  posted.\n"),
"user decision: the TMB internals note in the backlog")
cat("RELEASE FIX 2b DONE\n")

# readers closed the periodic-body minor after its confirmation (steps at
# irrational multiples, findings "Confirmation minor"); on the merged
# build the 10 pi row reads 0.8826 under frmtmb_mean_unsettled
# (dev/rel080-log/rr/readers-sc3.txt), so the backlog item says that
edit("dev/test-backlog.md",
paste0("- **MINOR. `fitted()` under a draw on a periodic body can drop the draw:**\n",
"  `y ~ sin(a)^2 + c` at `a = 0` with a draw sd of `10 pi` reads 0.509\n",
"  against 1.009, silent, because the steps 0.5, 1 and 1.7 sds land on zeros\n",
"  of the sine. Steps at irrational multiples would close it\n",
"  (dev/reviews/2026-11-06-readers.md, confirmation). Q1 of that review"),
paste0("- **MINOR. `fitted()` under a draw on a periodic body:** `y ~ sin(a)^2 + c`\n",
"  at `a = 0` with a draw sd of `10 pi` read 0.509 against 1.009, silent,\n",
"  on the reviewed build (dev/reviews/2026-11-06-readers.md, confirmation).\n",
"  The lane then stepped at irrational multiples (unreviewed); on 0.80.0\n",
"  the row reads 0.8826 under `frmtmb_mean_unsettled`. The rule cannot\n",
"  integrate `cos(20 pi z)`. Q1 of that review"),
"readers' later fix of the periodic-body minor")
cat("RELEASE FIX 2c DONE\n")

# nestquad closed its three review minors after the short check
# (findings "Fix round 2", unreviewed): the backlog keeps what remains
edit("dev/test-backlog.md",
paste0("  70 of 72 ones); the warning now gives the rule's own error. The memory\n",
"  fallback factor 42.9 over-predicts binomial and lognormal by 4.7 to 5.9\n",
"  times, so the default `quad_max_gb = 4` refuses them at 2000 to 4000\n",
"  rows where they need under 1 GB: calibrate them. `profile(fit)` with no\n",
"  `parm` stops on the missing argument before the nested refusal\n",
"  (dev/reviews/2026-11-06-nestquad.md, remaining items)."),
paste0("  70 of 72 ones). The lane then made the warning give the rule's own\n",
"  error (0.043 on the near-separable design, against 1.7e-2 measured),\n",
"  calibrated binomial, lognormal, von_mises and beta_binomial (predicted\n",
"  over measured memory 1.22 to 1.61 over 20 families) and refused\n",
"  `profile(fit)` with no `parm` by class, after the short check\n",
"  (dev/nestquad-findings.md, \"Fix round 2\"; not reviewed). What remains:\n",
"  the error itself (dev/reviews/2026-11-06-nestquad.md, remaining items)."),
"nestquad's later fix of its minors")
cat("RELEASE FIX 2d DONE\n")

# brief step 1: fitedge's replacement for the sentence of frm()'s
# `@param start` on how a nonlinear body starts (arstep rewrote the lines;
# dev/fitedge-findings.md, "Fix round", M2)
edit("R/fit.R",
paste0("#'   two. A nonlinear body, a mixture and `sigma ~ gp()` start at\n",
"#'   the best of a few length scales read from the objective. A `theta` here"),
paste0("#'   two. A field in a parameter of a nonlinear body starts at the\n",
"#'   marginal of the body made linear at its start, and the fit also\n",
"#'   descends from the best of a few length scales and keeps the better; a\n",
"#'   `start` of your own turns the first off. A mixture and `sigma ~ gp()`\n",
"#'   start at the best of a few length scales read from the objective. A\n",
"#'   `theta` here"),
"step 1, fitedge's sentence in frm()'s @param start")
cat("RELEASE FIX 2e DONE\n")
edit("R/fit.R",
paste0("#'   `theta` here\n",
"#'   replaces all of this, and a `beta` or `betad` replaces only itself.\n"),
paste0("#'   `theta` here replaces all of this, and a `beta` or `betad` replaces\n",
"#'   only itself.\n"),
"step 1, rewrap")

# bigx x nlstart, found on the merged build (dev/rel080-nlrows.R): when
# nl_basis_descent() kept its second descent, the tape carries a basis that
# whitens the body's Jacobian at the first fit's stop. The sign probe's
# patterns start at other signs, where that basis is wrong, and no pattern
# won: three of the review's fresh rows read 1.6 to 43.7 below the nlstart
# build, and 0 rows of 335 took the probe there. The runs step in the
# design's coordinates alone (the basis fit_assembled() starts from), and
# the fit's basis is put back after them.
edit("R/par-template.R",
"  saved <- obj_state_save(obj)
  best <- fit$opt
  for (s in pats) {",
"  saved <- obj_state_save(obj)
  # a basis from nl_basis_descent() whitens the body's Jacobian at the
  # first fit's stop; a pattern starts at other signs, where it does not
  # hold, and with it no pattern won on rows the probe exists for
  # (dev/rel080-nlrows.R). The runs step in the design's coordinates.
  ob_fit <- obj$env$.frmtmb_opt_basis
  if (isTRUE(ob_fit$nl)) {
    obj$env$.frmtmb_opt_basis <- opt_basis(frame, obj, ctl$optimizer)
  }
  best <- fit$opt
  for (s in pats) {",
"bigx x nlstart: the probe runs without the body's basis")
edit("R/par-template.R",
"  obj_settle(obj, fit$opt, fit$opt, saved)
  gain <- fit$opt$objective - best$objective",
"  obj$env$.frmtmb_opt_basis <- ob_fit
  obj_settle(obj, fit$opt, fit$opt, saved)
  gain <- fit$opt$objective - best$objective",
"bigx x nlstart: the fit's basis put back")
cat("RELEASE FIX 2f DONE\n")

# the test of the fix above; it fails on the build before it (the probe
# never wins once the body's basis is on the tape)
edit("tests/testthat/test-nl-start.R",
"  got <- caught(frm(sat_body, data = d))
  expect_true(has_class(got$conds, \"frmtmb_nl_multistart_unsettled\"))
  expect_true(has_class(got$conds, \"frmtmb_nl_multistart\"))
})",
"  got <- caught(frm(sat_body, data = d))
  expect_true(has_class(got$conds, \"frmtmb_nl_multistart_unsettled\"))
  expect_true(has_class(got$conds, \"frmtmb_nl_multistart\"))
})

test_that(\"the sign probe runs without a basis whitened at the first stop\", {
  # a Gompertz curve whose default fit stops at a lower maximum, where
  # nl_basis_descent() keeps a basis of the body's Jacobian; with that
  # basis the patterns' runs reached nothing (0.80.0 consolidation, 43.7
  # below the probe's point, review of lane nlstart's fresh row gomp 22)
  set.seed(22004)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  # the draws of the review's design, in its order
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  d$y <- 3 * exp(-2 * exp(-log(2) * d$xq)) + rnorm(n, 0, 0.5)
  fm <- bf(y ~ Asym * exp(-b2 * exp(-b3 * xq)), Asym ~ 1, b2 ~ 1, b3 ~ 1,
           nl = TRUE)
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             caught(frm(fm, data = d)))$value
  ob <- off$obj$env$.frmtmb_opt_basis
  # which point the first fit stops at depends on the optimizer's path;
  # the defect needs the body's basis on the tape
  skip_if_not(isTRUE(ob$nl), \"the first fit kept no basis of the body\")
  on <- caught(frm(fm, data = d))
  l_off <- as.numeric(logLik(off))
  l_on <- as.numeric(logLik(on$value))
  expect_gt(l_on - l_off, 1e-3 * abs(l_off))
  expect_true(has_class(on$conds, \"frmtmb_nl_multistart\"))
})",
"the test of the probe without the body's basis")
cat("RELEASE FIX 2g DONE\n")

# OpenBLAS 0.3.26 and 0.3.32 with LAPACK (r080-ob26, r080-ob32), two tests
# that read one optimizer path:
# test-se-reparam.R: the bounded degree-6 fit stops at code 9 (evaluation
# limit) at the maximum under OpenBLAS, and says so with the Newton step's
# gain (6.81e-08) within grad_tol; the test already checks the maximum to
# 1e-7, so that stop is a valid outcome
edit("tests/testthat/test-se-reparam.R",
"    frm(fo, data = d, prior = set_prior(\"\", class = \"b\", coef = \"x6\",
                                        ub = ub)),
    \"x6\", require = \"x6\")",
"    frm(fo, data = d, prior = set_prior(\"\", class = \"b\", coef = \"x6\",
                                        ub = ub)),
    # under OpenBLAS the fit can stop at its evaluation limit at the
    # maximum, and says so; the log-likelihood below is the check
    c(\"x6\", \"did not report convergence\"), require = \"x6\")",
"OpenBLAS: the bounded fit may stop at code 9 at the maximum")
# test-nl-start.R: the estimate of b2 agrees to 1.7e-4 relative under
# OpenBLAS at a maximum the log-likelihood check holds to 1e-6; the
# estimate is compared in units of its standard error
edit("tests/testthat/test-nl-start.R",
"  # one model with b2 rescaled
  expect_equal(at6$value$estimates$beta[[\"b2_(Intercept)\"]] * 1e-6,
               at1$estimates$beta[[\"b2_(Intercept)\"]], tolerance = 1e-4)",
"  # one model with b2 rescaled, to a small part of b2's standard error
  # (under OpenBLAS the two estimates part by 1.7e-4 relative, at a
  # maximum the log-likelihood check above holds)
  fe <- fixef(at1)
  se_b2 <- fe[grep(\"^b2_\", rownames(fe)), \"Est.Error\"]
  expect_lt(abs(at6$value$estimates$beta[[\"b2_(Intercept)\"]] * 1e-6 -
                  at1$estimates$beta[[\"b2_(Intercept)\"]]), 0.01 * se_b2)",
"OpenBLAS: b2 compared in SE units")
cat("RELEASE FIX 2h DONE\n")
