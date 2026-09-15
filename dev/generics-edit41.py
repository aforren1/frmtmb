# The review's sentence about the guard: the bare-UseMethod test stops
# work coming BACK into the generic, and cannot see work that LEFT the
# generic and missed a method. This is the other half: every
# hypothesis() method in core and in frmtmb.sample arms the note.
P = ("C:/Users/adf44/source/r/frmtmb-wt-generics/extensions/frmtmb.sample/"
     "tests/testthat/test-draws-methods.R")
s = open(P, encoding="utf-8", newline="").read()
marker = 'test_that("every hypothesis() method arms the note itself'
assert marker not in s, "already appended"
block = '''
test_that("every hypothesis() method arms the note itself", {
  # The structural half of the test above, and the reason it exists.
  #
  # Core guards that no borrowed generic carries work in its body
  # (tests/testthat/test-generic-collision.R), because the exported
  # generic may be brms's and then that work does not run. That guard
  # stops work coming BACK into a generic. It cannot see work that LEFT
  # the generic and did not reach every method, which is exactly what
  # happened here: the arming moved into core's two methods and this
  # package's draws method was not one of them. So this asserts where
  # the work WENT, over every hypothesis() method both packages define.
  meths <- c(
    ls(asNamespace("frmtmb"), all.names = TRUE, pattern = "^hypothesis[.]"),
    ls(asNamespace("frmtmb.sample"), all.names = TRUE,
       pattern = "^hypothesis[.]"))
  # the guard is only a guard if it found the methods it is about
  expect_true(all(c("hypothesis.frmtmb_fit", "hypothesis.frmtmb_multiple",
                    "hypothesis.frmtmb_draws") %in% meths))
  arms <- vapply(meths, function(m) {
    ns <- if (exists(m, envir = asNamespace("frmtmb.sample"),
                     inherits = FALSE)) "frmtmb.sample" else "frmtmb"
    f <- get(m, envir = asNamespace(ns), inherits = FALSE)
    b <- paste(deparse(body(f)), collapse = " ")
    grepl("hyp_shadow_arm()", b, fixed = TRUE) &&
      grepl("hyp_shadow_disarm(", b, fixed = TRUE)
  }, NA)
  expect_equal(names(arms)[!arms], character())
})
'''
open(P, "w", encoding="utf-8", newline="\n").write(s + block)
print("appended")
