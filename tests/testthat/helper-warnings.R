# Evaluate `expr`, muffling only the warnings whose message contains one
# of the strings in `allowed`. Every other warning still escapes to
# testthat, which counts and prints it.
#
# Why not expect_warning(): in testthat 3 it absorbs the FIRST warning
# it matches and lets the rest escape, so a fit that warns ten times on
# purpose left nine warnings in the suite's WARN count, beside the ones
# nobody expected (dev/warnleak-scan.R lists them). Why not
# suppressWarnings(): it would also hide a warning nobody expected.
#
# `require` names messages that MUST occur, each at least once, for a
# test whose point is that the warning is raised. The expression is
# evaluated lazily in the caller's frame, so `allow_warnings(fit <-
# frm(...), ...)` assigns `fit` there.
# The warning of class frmtmb_interop_se_part: marginaleffects on a fit with a
# variance no parameter of get_coef() holds (a me() draw, a latent residual, a
# fit that integrates its group-level effects out). A test whose point is
# something else allows it.
mfx_se_part <- "standard errors from the parameters get_coef() holds"

allow_warnings <- function(expr, allowed, require = character()) {
  seen <- character()
  val <- withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (any(vapply(allowed, function(a) grepl(a, m, fixed = TRUE), NA))) {
      seen <<- c(seen, m)
      invokeRestart("muffleWarning")
    }
  })
  for (r in require) {
    testthat::expect_true(any(grepl(r, seen, fixed = TRUE)),
                          info = paste0("expected a warning containing '",
                                        r, "'"))
  }
  invisible(val)
}

# Evaluate `expr`, muffling the message of class frmtmb_converged_at_max
# (nlminb's "false convergence (8)" at a measured maximum) and nothing
# else, for a test whose fit may stop there on some rounding and whose
# point is something else.
allow_at_max <- function(expr) {
  withCallingHandlers(expr, frmtmb_converged_at_max = function(m) {
    invokeRestart("muffleMessage")
  })
}

# Evaluate `expr`, muffling the boundary (singular) fit message of the
# standard-error check (class frmtmb_boundary_fit) and nothing else, for
# a test whose data happen to put a variance component at zero and whose
# point is some other condition. `require = TRUE` asserts the message
# was given.
allow_boundary <- function(expr, require = FALSE) {
  seen <- 0L
  val <- withCallingHandlers(expr, frmtmb_boundary_fit = function(m) {
    seen <<- seen + 1L
    invokeRestart("muffleMessage")
  })
  if (require) {
    testthat::expect_gt(seen, 0L)
  }
  invisible(val)
}
