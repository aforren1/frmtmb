
# prior_summary() of draws lists the rows the user wrote, both spellings,
# as frmtmb.sample 0.24.0 and brms do, while the sampler uses the merged
# prior. Release fix R1 printed the one merged class row in their place:
# `normal(5, 0.1)` on the class with the class row's bound, a row nobody
# wrote, and the user's `normal(10, 1)` was gone
# (dev/reviews/2026-10-29-release.md, finding 2).
test_that("prior_summary() keeps both spellings while the sampler merges them", {
  skip_sampler()
  d <- pni_data()
  pA <- set_prior("normal(5, 0.1)", class = "Intercept", nlpar = "a")
  pB2 <- set_prior("normal(10, 1)", class = "b", nlpar = "a", lb = 2)
  fit <- pni_quiet(frm(pni_form(), data = d, prior = pA))
  ds <- pni_quiet(frm_sample(fit, prior = pB2, chains = 1, iter = 100,
                             seed = 3, refresh = 0))
  tab <- as.data.frame(prior_summary(ds))
  ra <- tab[tab$nlpar == "a", c("prior", "class", "coef", "lb")]
  # the two rows as written: the coefficient row and the class row
  expect_identical(nrow(ra), 2L)
  expect_identical(sort(ra$prior), c("normal(10, 1)", "normal(5, 0.1)"))
  expect_identical(ra$coef[ra$prior == "normal(5, 0.1)"], "Intercept")
  expect_identical(ra$coef[ra$prior == "normal(10, 1)"], "")
  expect_identical(as.numeric(ra$lb[ra$prior == "normal(10, 1)"]), 2)
  # the sampler's prior is the merged one: the intercept spelling's
  # density with the class row's bound
  ri <- pni_ri(fit, pA + pB2)
  a <- pni_a(ri, fit)
  expect_identical(a$dist[[1L]]$location, 5)
  expect_identical(a$lower, 2)
  # and the printed rows, passed back, resolve to that same prior
  expect_identical(pni_a(pni_ri(fit, prior_summary(ds)), fit), a)
})
