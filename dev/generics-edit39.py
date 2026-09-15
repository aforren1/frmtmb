# BLOCKER E, step 1: the frmtmb.sample test, appended BEFORE the fix so
# it can be seen failing against the current tree.
P = ("C:/Users/adf44/source/r/frmtmb-wt-generics/extensions/frmtmb.sample/"
     "tests/testthat/test-draws-methods.R")
s = open(P, encoding="utf-8", newline="").read()
marker = 'test_that("hypothesis() on draws gives the reserved-name note'
assert marker not in s, "already appended"
block = '''
## ---- hypothesis() naming notes --------------------------------------

test_that("hypothesis() on draws gives the reserved-name note, once", {
  # A covariate named `sigma` shadows the residual SD, and hypothesis()
  # says which one it read. The note is ARMED per user-level call and
  # emitted deep inside hyp_env_vals(), so the method the user reached
  # has to arm it.
  #
  # Core used to arm it in its GENERIC. Core's exported `hypothesis`
  # now resolves to brms's generic whenever brms is loaded, so the
  # arming moved into core's own methods, and this method, which had
  # relied on the generic, lost the note in EVERY session, brms or not.
  # Measured before the fix: 1 note on the frmtmb_fit, 0 on its draws.
  #
  # fake_draws() rather than the sampler: the note is emitted while the
  # hypothesis is parsed against the fit, before any draw is read, so
  # zero draws exercise exactly the path in question and this test
  # needs no Stan build to run.
  set.seed(3)
  n <- 200
  dd <- data.frame(sigma = stats::rnorm(n),
                   g = factor(rep(1:10, length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.7 * dd$sigma +
                         stats::rnorm(10, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ sigma + (1 | g)), family = gaussian(), data = dd)
  pat <- "reads 'sigma' as the coefficient"

  # the control, on the fit, so a missing note on draws cannot be a
  # construction that never shadows anything
  on_fit <- capture_messages(hypothesis(fit, "sigma = 0"))
  expect_equal(sum(grepl(pat, on_fit, fixed = TRUE)), 1L)

  ds <- fake_draws(fit)
  on_draws <- capture_messages(
    suppressWarnings(hypothesis(ds, "sigma = 0")))
  expect_equal(sum(grepl(pat, on_draws, fixed = TRUE)), 1L)
})
'''
if not s.endswith("\n"):
    s += "\n"
open(P, "w", encoding="utf-8", newline="\n").write(s + block)
print("appended")
