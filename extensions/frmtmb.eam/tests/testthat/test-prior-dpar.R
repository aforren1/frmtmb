# An accumulator model's location is one drift rate per accumulator
# (v1, v2, ...), and frmtmb addresses each by dpar, as brms does for the
# several-location families it has. Through frmtmb 0.62.0 a class "b" or
# "Intercept" prior with no dpar reached every drift at once, and
# default_prior() listed one row with an empty dpar
# (dev/mvprior-log/extfam-base.txt, -lane.txt in the frmtmb tree, script
# dev/mvprior-extfam.R).

pd_msg <- function(expr) {
  e <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (inherits(e, "frmtmb_error")) conditionMessage(e) else ""
}

pd_data <- function(sim) {
  set.seed(1)
  d <- sim(300, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.45, ndt = 0.2)
  d$x <- stats::rnorm(nrow(d))
  d
}

test_that("lba() and rdm(): b and Intercept with no dpar are refused", {
  for (case in list(list(lba(3), lba_simulate), list(rdm(3), rdm_simulate))) {
    d <- pd_data(case[[2]])
    for (p in list(set_prior("normal(0, 1)", class = "b"),
                   set_prior("normal(0, 1)", class = "Intercept"))) {
      expect_match(pd_msg(validate_prior(p, bf(rt | vint(choice) ~ x),
                                         data = d, family = case[[1]])),
                   "several distributional parameters, v1, v2, v3")
    }
  }
})

test_that("default_prior() lists each drift by dpar", {
  tab <- as.data.frame(default_prior(bf(rt | vint(choice) ~ x),
                                     data = pd_data(lba_simulate),
                                     family = lba(3)))
  expect_true(all(c("v1", "v2", "v3") %in% tab$dpar[tab$class == "b"]))
  expect_false(any(tab$class %in% c("b", "Intercept") & !nzchar(tab$dpar)))
})

# Passes on 0.62.0 too: the dpar spelling reaches one drift.
test_that("a prior with dpar reaches that drift only", {
  des <- frmtmb:::prior_design(bf(rt | vint(choice) ~ x),
                               pd_data(rdm_simulate), rdm(3), list())
  r <- frmtmb:::resolve_priorlist(des, set_prior("normal(0, 1)", class = "b",
                                                 dpar = "v2"))
  pt <- des$frame[["par_template"]]
  nm <- unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  }))
  expect_identical(nm, "v2_x")
})
