## The scale guard. Every default in this package is written for a
## response time in seconds, and nothing in any of the five likelihoods
## refuses milliseconds: the fit converges and reports a boundary
## separation three orders of magnitude out. These tests pin the
## warning that says so, on all five families, and pin the band it
## fires in.

# One deterministic set of response times on the scale the package
# expects. Reused so that the seconds arm and the millisecond arm of
# every check differ in nothing but the factor of 1000.
sec_rt <- c(0.281, 0.334, 0.402, 0.455, 0.517, 0.588, 0.664, 0.751,
            0.849, 0.962, 1.093, 1.248, 1.435, 1.667, 1.964, 2.371)
ms_rt <- sec_rt * 1000

# Collect only the units warning, muffling everything else, so that a
# family raising an unrelated warning cannot make one of these pass or
# fail for the wrong reason.
units_warnings <- function(expr) {
  got <- character(0)
  withCallingHandlers(
    force(expr),
    frmtmb_eam_units_warning = function(w) {
      got <<- c(got, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  got
}

# The five families and the addition terms each one's valid_y reads,
# built at whatever scale the response is on so that a family with its
# own scale-bearing argument (the go/no-go deadline) stays consistent.
eam_cases <- function(rt) {
  n <- length(rt)
  up <- rep(0:1, length.out = n)
  list(
    wiener = list(fam = wiener(),
                  aterms = list(dec = up)),
    lba = list(fam = lba(2),
               aterms = list(vint1 = up + 1)),
    rdm = list(fam = rdm(2),
               aterms = list(vint1 = up + 1)),
    gddm = list(fam = gddm(),
                aterms = list(dec = up, vint1 = rep(1L, n))),
    wiener_gng = list(fam = wiener_gng(deadline = 4 * max(rt)),
                      aterms = list(dec = up)))
}

test_that("every family warns on a millisecond response", {
  cases <- eam_cases(ms_rt)
  for (nm in names(cases)) {
    got <- units_warnings(
      cases[[nm]]$fam[["valid_y"]](ms_rt, cases[[nm]]$aterms))
    expect_length(got, 1)
    # the family the user wrote, not the helper that raised it
    expect_match(got, nm, fixed = TRUE)
    expect_match(got, "millisecond", ignore.case = TRUE)
    # what to divide by, which is the only actionable part
    expect_match(got, "1000", fixed = TRUE)
  }
})

test_that("no family warns on a seconds response", {
  cases <- eam_cases(sec_rt)
  for (nm in names(cases)) {
    got <- units_warnings(
      cases[[nm]]$fam[["valid_y"]](sec_rt, cases[[nm]]$aterms))
    expect_length(got, 0)
  }
})

test_that("the band is the documented one, and open at the bottom", {
  ceil <- ddm_seconds_ceiling
  up <- rep(0:1, length.out = 4L)
  # a design whose fastest response sits exactly at the ceiling is
  # taken at its word: the guard reads "above", not "at or above"
  at <- ceil * c(1, 1.5, 2, 3)
  expect_length(units_warnings(wiener()[["valid_y"]](at, list(dec = up))),
                0)
  # and one whose fastest response is the smallest step above it warns
  just_over <- at
  just_over[1L] <- ceil * (1 + 8 * .Machine$double.eps)
  expect_gt(just_over[1L], ceil)
  expect_length(
    units_warnings(wiener()[["valid_y"]](just_over, list(dec = up))), 1)
})

test_that("the warning carries its own class", {
  # a deliberation task with no trial under the ceiling is the false
  # alarm this guard is allowed to make, and such a design has to be
  # able to silence THIS condition without also silencing the
  # convergence warnings beside it.
  cond <- NULL
  withCallingHandlers(
    wiener()[["valid_y"]](ms_rt, list(dec = rep(0:1, 8))),
    warning = function(w) {
      cond <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(cond, "frmtmb_eam_units_warning")
  expect_s3_class(cond, "warning")
})

test_that("the warning reaches a user through frm()", {
  # reachability, not a fit: valid_y runs at frame assembly, so one
  # optimizer step is enough to prove the path a user actually takes.
  set.seed(414)
  n <- 40L
  dat <- data.frame(rt = rep(ms_rt, length.out = n),
                    upper = rep(0:1, length.out = n))
  got <- units_warnings(
    try(frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(),
            data = dat,
            control = frmtmb_control(optCtrl = list(iter.max = 1L,
                                                    eval.max = 2L))),
        silent = TRUE))
  expect_length(got, 1)
  expect_match(got, "wiener", fixed = TRUE)
})
