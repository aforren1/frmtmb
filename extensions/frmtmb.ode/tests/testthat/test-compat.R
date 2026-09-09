# What this package tells frmtmb's compatibility registry about itself.
#
# The registration happens in .onLoad(), so by the time this file runs
# it has already happened; what is asserted here is that it took.
#
# The rule refuses `frm_ode() x frm_sample`. That pair only exists when
# frmtmb.sample is loaded, so the rule is declared through
# `expects = "frm_sample"` and lies dormant otherwise. `frm_compat()`
# therefore cannot be asked about the pair here, while
# `frm_compat_rules()` reports the rule itself either way; the end-to-
# end check that frm_sample() reads the row lives in that package's
# tests/testthat/test-compat-preflight.R.

test_that("frm_ode() is in the compatibility vocabulary", {
  ft <- frmtmb::frm_compat_features()
  expect_true("frm_ode()" %in% ft$name)
  expect_equal(ft$kind[ft$name == "frm_ode()"], "special")
  # the key is what a caller matches against a call in a formula
  expect_equal(ft$key[ft$name == "frm_ode()"], "frm_ode")
})

test_that("sampling an frm_ode() model is registered as refused", {
  rl <- frmtmb::frm_compat_rules()
  row <- rl[rl$feature_a == "frm_ode()" & rl$feature_b == "frm_sample", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$status, "refused")
  # the note is what frm_sample() repeats back, so it has to say what
  # to do instead rather than only what is wrong
  expect_match(row$note, "maximum likelihood")
  expect_match(row$note, "dev/upstream")
})

test_that("the registration survives being run again, which is what
           pkgload::load_all() does to a namespace", {
  # frmtmb_register_compat(features =) is a no-op for a name it already
  # carries, so a second run must not error. The RULES are appended
  # rather than deduplicated, which is core's behavior and not this
  # package's: frmtmb.sample's six rules double the same way. Identical
  # rules tie on specificity and agree on the status, so the resolved
  # answer is unchanged, which is what this asserts.
  #
  # This block LEAVES a duplicate rule in the registry for the rest of
  # the R process, and there is no way to withdraw one. Under
  # `R CMD check` testthat runs this package's files in one process,
  # and this file sorts before every test-ode*.R, none of which reads
  # the registry, so nothing downstream is affected. Anyone adding a
  # registry-reading file to this package should either put it before
  # this one or expect the extra rule.
  before <- nrow(frmtmb::frm_compat_rules())
  expect_silent(frmtmb.ode:::.onLoad(NULL, "frmtmb.ode"))
  rl <- frmtmb::frm_compat_rules()
  row <- rl[rl$feature_a == "frm_ode()" & rl$feature_b == "frm_sample", ]
  expect_gte(nrow(row), 1L)
  expect_setequal(row$status, "refused")
  expect_gte(nrow(rl), before)
})
