# The default-prior announcement for a simplex a family declares
# (frmtmb_family(simplex_pars =)). It used to name cox()'s baseline by
# family name, so a third-party family's simplex took brms's
# dirichlet(1) on the sampler without the announcement saying so. The
# family here is frmtmb's test-simplex-pars.R one: a location mixture
# over three fixed centers with simplex weights.

tsimplex_family <- function(centers = c(-3, 0, 3)) {
  frmtmb_family(
    "tsimplex",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      w <- frmtmb::mo_simplex(extra$tw_raw, aterms[["tw_chart"]])
      cen <- aterms[["tw_centers"]]
      f <- 0
      for (j in seq_along(cen)) {
        f <- f + w[j] * exp(-0.5 * (y - dpars[["mu"]] - cen[j])^2)
      }
      log(f) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) 0),
    extra_pars = function(y, aterms) {
      list(tw_raw = rep(0, length(aterms[["tw_centers"]]) - 1L))
    },
    aterm_data = function(y, aterms) list(tw_centers = centers),
    simplex_pars = c(tw = "tw_raw"))
}

tsimplex_data <- function(seed = 5, n = 200) {
  set.seed(seed)
  k <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  data.frame(y = 0.5 + c(-3, 0)[k] + stats::rnorm(n))
}

test_that("a declared simplex's default prior is announced by its name", {
  fit <- frm(y ~ 1, family = tsimplex_family(), data = tsimplex_data())
  n <- default_prior_notes(fit)
  hit <- grepl("tw: dirichlet(1), brms's default", n, fixed = TRUE)
  expect_identical(sum(hit), 1L)
  expect_match(n[hit], "the tsimplex() family declares", fixed = TRUE)
  # the guard absent: a user's row on the class takes the density, and
  # the announcement steps aside
  pr <- set_prior("dirichlet(20)", class = "tw")
  expect_false(any(grepl("tw: dirichlet(1)", default_prior_notes(fit, pr),
                         fixed = TRUE)))
  # and a model with no declared simplex announces none
  g <- frm(y ~ 1, family = gaussian(), data = tsimplex_data())
  expect_false(any(grepl("family declares", default_prior_notes(g),
                         fixed = TRUE)))
})
