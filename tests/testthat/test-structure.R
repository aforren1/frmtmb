# The structured-family protocol's constructor and the core accessors
# that read it. The families that RIDE the protocol are covered by
# test-hmm.R and the mixture tests; what is asserted here is the
# contract a family written outside this package must meet.

ll_ok <- function(y, dpars, aterms, weights, block, extra) sum(y)

test_that("every slot defaults to the rowwise behavior", {
  expect_error(frmtmb_structure(loglik = "not a function"), "loglik")
  st <- frmtmb_structure(loglik = ll_ok)
  expect_s3_class(st, "frmtmb_structure")
  expect_null(st[["frame_vars"]])
  expect_null(st[["fitted_mean"]])
  expect_false(st[["keep_na"]])
})

test_that("a structure with no loglik declares capabilities only", {
  # the lca() shape: a likelihood that DOES factorize per row, carried
  # by the family's own lpdf, plus refusals that belong to the family
  st <- frmtmb_structure(supports = list(osa = FALSE))
  expect_s3_class(st, "frmtmb_structure")
  expect_null(st[["loglik"]])
  expect_false(structure_allows(st, "osa"))
})

test_that("every capability starts refused", {
  st <- frmtmb_structure(loglik = ll_ok)
  expect_false(any(st[["supports"]]))
  expect_length(st[["supports"]], length(frmtmb_structure_flags))
  # a family with no structure at all is rowwise and refuses nothing
  expect_true(structure_allows(NULL, "osa"))
  expect_false(structure_allows(st, "osa"))
})

test_that("unit is one noun phrase, or the generic one", {
  expect_error(frmtmb_structure(loglik = ll_ok, unit = c("a", "b")),
               "noun phrase")
  expect_error(frmtmb_structure(loglik = ll_ok, unit = ""), "noun phrase")
  expect_identical(
    structure_unit(frmtmb_structure(loglik = ll_ok, unit = "a sequence")),
    "a sequence")
  expect_true(nzchar(structure_unit(frmtmb_structure(loglik = ll_ok))))
})

test_that("the constructor validates slot types", {
  for (arg in c("frame_vars", "check_spec", "frame_block", "check_frame",
                "check_fit", "fitted_mean", "fitted_var", "latent_probs",
                "sim_ctx")) {
    a <- list(loglik = ll_ok)
    a[[arg]] <- 1
    expect_error(do.call(frmtmb_structure, a), "must be a function")
  }
  expect_error(frmtmb_structure(loglik = ll_ok, keep_na = NA), "TRUE or FALSE")
  expect_error(frmtmb_structure(loglik = ll_ok, keep_na = "yes"),
               "TRUE or FALSE")
})

test_that("supports names come from a closed vocabulary", {
  expect_error(frmtmb_structure(loglik = ll_ok, supports = list(nope = TRUE)),
               "unknown capability flag")
  expect_error(frmtmb_structure(loglik = ll_ok, supports = list(TRUE)),
               "must name every flag")
  expect_error(frmtmb_structure(loglik = ll_ok, supports = list(osa = "yes")),
               "must be TRUE or FALSE")
  expect_error(
    frmtmb_structure(loglik = ll_ok,
                     supports = list(osa = TRUE, osa = FALSE)),
    "names a flag twice")
  st <- frmtmb_structure(loglik = ll_ok, supports = c(osa = TRUE))
  expect_true(st[["supports"]][["osa"]])
  expect_false(st[["supports"]][["reml"]])
})

test_that("a refusal must explain a flag that is actually refused", {
  expect_error(frmtmb_structure(loglik = ll_ok, refusals = list(zz = "x")),
               "explains no known capability flag")
  expect_error(
    frmtmb_structure(loglik = ll_ok, supports = list(osa = TRUE),
                     refusals = list(osa = "x")),
    "could never be shown")
  expect_error(frmtmb_structure(loglik = ll_ok, refusals = list(osa = 1)),
               "one non-empty string")
  expect_error(frmtmb_structure(loglik = ll_ok, refusals = list("x")),
               "must name the flag")
})

test_that("a loglik that drops the weights the core passes is named", {
  expect_warning(
    frmtmb_structure(loglik = function(y, dpars, aterms, block, extra) 1),
    "weights")
  expect_silent(frmtmb_structure(loglik = ll_ok))
  # a dots-taking loglik has not dropped anything
  expect_silent(frmtmb_structure(loglik = function(y, ...) 1))
})

test_that("a refusal is raised in the family's own words", {
  st <- frmtmb_structure(
    loglik = ll_ok,
    refusals = list(osa = "no one-step residual here",
                    re_form = "no re.form here",
                    re_form.simulate = "not when simulating"))
  expect_error(structure_gate(st, "osa", "GENERIC"),
               "no one-step residual here", fixed = TRUE)
  # a context-specific message wins, and falls back to the bare flag
  expect_error(structure_gate(st, "re_form", "GENERIC",
                              context = "simulate"),
               "not when simulating", fixed = TRUE)
  expect_error(structure_gate(st, "re_form", "GENERIC"),
               "no re.form here", fixed = TRUE)
  # a flag refused without a sentence of its own gets the generic one
  expect_error(structure_gate(st, "deviance", "GENERIC"), "GENERIC",
               fixed = TRUE)
  # nothing is raised for a supported flag or for a rowwise family
  expect_silent(structure_gate(NULL, "osa", "GENERIC"))
})

test_that("frmtmb_family() takes a structure and nothing else", {
  st <- frmtmb_structure(loglik = ll_ok)
  fam <- frmtmb_family("demo", "mu", list(mu = "identity"),
                       lpdf = function(y, dpars, aterms) 0, structure = st)
  expect_identical(fam_structure(fam), st)
  expect_null(fam_structure(as_frmtmb_family(stats::gaussian())))
  expect_error(
    frmtmb_family("demo", "mu", list(mu = "identity"),
                  lpdf = function(y, dpars, aterms) 0,
                  structure = list(loglik = ll_ok)),
    "frmtmb_structure")
})

test_that("a group-level mixture carries the protocol", {
  fg <- mixture(gaussian(), gaussian(), groups = ~g)
  st <- fam_structure(fg)
  expect_s3_class(st, "frmtmb_structure")
  # the four flags that follow from the per-row mean being rowwise
  for (flag in c("conditional_effects", "newdata_response",
                 "se_fit_response", "re_form")) {
    expect_true(st[["supports"]][[flag]], label = flag)
  }
  # the group branch registers no observation vector, so there is
  # nothing for one-step-ahead residuals to step through
  expect_false(st[["supports"]][["osa"]])
  expect_true(is.function(st[["loglik"]]))
})

test_that("a rowwise mixture-type family declares capabilities only", {
  # mixture(), mixture_mvn() and lca() all have a likelihood that DOES
  # factorize per row, so their structures carry no loglik; what they
  # carry is the multimodality refusal that used to be one gate in
  # fit.R naming all three
  for (fam in list(mixture(gaussian(), gaussian()),
                   mixture_mvn(K = 2, D = 2), lca(K = 2))) {
    st <- fam_structure(fam)
    expect_s3_class(st, "frmtmb_structure")
    expect_null(st[["loglik"]])
    expect_false(st[["supports"]][["reml"]])
    expect_false(st[["supports"]][["profile"]])
    # everything else stays exactly as available as it was
    expect_true(st[["supports"]][["quadrature"]])
    expect_true(st[["supports"]][["cluster_robust"]])
    expect_true(is.function(st[["latent_probs"]]))
  }
  # only the lca refuses a one-step-ahead residual
  expect_false(fam_structure(lca(K = 2))[["supports"]][["osa"]])
  expect_true(
    fam_structure(mixture(gaussian(), gaussian()))[["supports"]][["osa"]])
})

test_that("an hmm() carries the protocol with every capability refused", {
  st <- fam_structure(hmm(K = 2, gaussian(), time = t, group = id))
  expect_s3_class(st, "frmtmb_structure")
  expect_false(any(st[["supports"]]))
  expect_true(st[["keep_na"]])
  for (slot in c("frame_vars", "check_spec", "frame_block", "loglik",
                 "fitted_mean", "fitted_var", "latent_probs", "sim_ctx")) {
    expect_true(is.function(st[[slot]]), label = slot)
  }
})

test_that("the frame carries one block per structured response", {
  set.seed(4)
  ng <- 12L
  m <- 5L
  g <- rep(seq_len(ng), each = m)
  cls <- stats::rbinom(ng, 1L, 0.5)
  dd <- data.frame(y = stats::rnorm(ng * m, c(-2, 2)[cls + 1L][g], 0.5),
                   g = factor(g))
  fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
             data = dd)
  blk <- fit$frame[["blocks"]][["y"]]
  expect_type(blk, "list")
  expect_identical(length(blk[["first"]]), ng)
  expect_identical(blk[["levels"]], levels(dd$g))
  # the block is data, so a refit rebuilds an identical one
  expect_null(fit$frame[["mix_g"]])

  # a rowwise fit carries no block at all
  f2 <- frm(bf(y ~ 1) + gaussian(), data = dd)
  expect_length(f2$frame[["blocks"]], 0L)
})
