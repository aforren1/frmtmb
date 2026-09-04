# Parallel chains reach draws on every platform. On Windows the chains
# run on PSOCK workers: the tape's external pointer dies in
# serialization, tmbstan retapes on the worker from the objective
# closure, and the closure must therefore be self-contained. This test
# is the structural regression gate for that property; it is what the
# old sequential-fallback guard wrongly claimed impossible.

test_that("frm_sample runs parallel chains, formula route with ncp", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  # a load_all() dev namespace cannot be rebuilt on a worker; only the
  # installed package serializes
  skip_if(exists(".__DEVTOOLS__", asNamespace("frmtmb")),
          "parallel workers need the installed package")
  dd <- data.frame(
    x = rnorm(80),
    g = factor(rep(1:8, 10))
  )
  set.seed(4)
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  # suppressWarnings: 150 post-warmup draws trip rstan's ESS and R-hat
  # advice, and this test is structural (did both workers deliver a
  # chain), not a mixing gate
  expr <- quote(suppressWarnings(
    frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
               data = dd, chains = 2, cores = 2, iter = 300,
               warmup = 150, seed = 7, refresh = 0)))
  ds <- NULL
  if (.Platform$OS.type == "windows") {
    expect_message(ds <- eval(expr), "parallel chains on Windows")
  } else {
    ds <- eval(expr)
  }
  # 2 chains x 150 post-warmup draws: a chain that died on its worker
  # would halve this
  expect_identical(posterior::ndraws(as_draws(ds)), 300L)
})
