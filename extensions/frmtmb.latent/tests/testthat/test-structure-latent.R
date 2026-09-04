# The two rows this package owns in frmtmb's structured-family table.
#
# frmtmb's own test-structure.R asserts the constructor's contract and
# the mixture rows of that table. These two blocks came from it in
# protocol step 10, when hmm() and lca() moved out: the assertions are
# unchanged, and the only edit is that a structure is read as
# `fam[["structure"]]` rather than through frmtmb's internal
# fam_structure(), which an extension cannot see and does not need.

test_that("a rowwise mixture-type family declares capabilities only", {
  # lca() has a likelihood that DOES factorize per row, so its
  # structure carries no loglik; what it carries is the multimodality
  # refusal that used to be one gate in frmtmb's fit.R naming every
  # mixture-type family.
  st <- lca(K = 2)[["structure"]]
  expect_s3_class(st, "frmtmb_structure")
  expect_null(st[["loglik"]])
  expect_false(st[["supports"]][["reml"]])
  expect_false(st[["supports"]][["profile"]])
  # everything else stays exactly as available as it was
  expect_true(st[["supports"]][["quadrature"]])
  expect_true(st[["supports"]][["cluster_robust"]])
  expect_true(is.function(st[["latent_probs"]]))
  # only the lca refuses a one-step-ahead residual
  expect_false(st[["supports"]][["osa"]])
})

test_that("an hmm() carries the protocol with every capability refused", {
  st <- hmm(K = 2, gaussian(), time = t, group = id)[["structure"]]
  expect_s3_class(st, "frmtmb_structure")
  expect_false(any(st[["supports"]]))
  expect_true(st[["keep_na"]])
  for (slot in c("frame_vars", "check_spec", "frame_block", "loglik",
                 "fitted_mean", "fitted_var", "latent_probs", "sim_ctx")) {
    expect_true(is.function(st[[slot]]), label = slot)
  }
})

test_that("both families reach frmtmb only through its exported seams", {
  # the out-of-tree compilation check as a test rather than a one-off:
  # every frmtmb symbol these two files call by name must be an export.
  # dev/out-of-tree-inventory.md is the written form of this scan, and
  # this is what fails if a later edit reaches for an internal again.
  skip_if_not_installed("codetools")
  ns <- asNamespace("frmtmb.latent")
  own <- ls(ns, all.names = TRUE)
  fns <- Filter(function(o) is.function(get(o, envir = ns)), own)
  used <- unique(unlist(lapply(fns, function(o)
    codetools::findGlobals(get(o, envir = ns), merge = TRUE))))
  core <- ls(asNamespace("frmtmb"), all.names = TRUE)
  reached <- setdiff(intersect(used, core), own)
  expect_setequal(setdiff(reached, getNamespaceExports("frmtmb")),
                  character(0))
})
