# A draws object with `stanfit = NULL` is what a test builds when it
# wants chosen parameter vectors without a sampler. The chain count was
# read as `x$stanfit@sim$chains %||% 1L`, and `@` on NULL is an error
# that never reaches `%||%`, so nchains(), summary(), VarCorr(), the
# as_draws_*() converters and every other reader of the chain count died
# on such an object (dev/sampfix-log/03-ref.txt). They count it as one
# chain now; the two readers of the sampler's own diagnostics refuse it
# by name.

nosf_case <- function(n = 40L) {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x + stats::rnorm(6, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
             dry_run = "objective")
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(stats::rnorm(n * length(lab), 0, 0.1) + 0.5, n,
              dimnames = list(NULL, lab))
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

test_that("a draws object without a stanfit counts as one chain", {
  ds <- nosf_case()
  expect_identical(nchains(ds), 1L)
  expect_identical(niterations(ds), 40L)
  expect_identical(ndraws(ds), 40L)
  expect_identical(dim(as_draws_array(ds))[1:2], c(40L, 1L))
  expect_identical(dim(as.array(ds))[1:2], c(40L, 1L))
  expect_s3_class(as_draws_df(ds), "draws_df")
  expect_s3_class(as_draws_rvars(ds), "draws_rvars")
  expect_identical(nrow(as.matrix(ds)), 40L)
  expect_identical(nrow(as.data.frame(ds)), 40L)
  expect_true(is.matrix(summary(ds)))
  expect_true(is.matrix(posterior_summary(ds)))
  expect_true(is.matrix(posterior_interval(ds)))
  expect_true(is.matrix(fixef(ds)))
  vc <- VarCorr(ds)
  expect_named(vc, c("g", "residual__"))
  expect_output(print(vc))
  expect_identical(attr(VarCorr(ds, summary = FALSE)$g$sd, "nchains"), 1L)
  expect_named(ranef(ds), "g")
  expect_true(is.numeric(rhat(ds)))
  skip_if_not_installed("coda")
  expect_s3_class(as.mcmc(ds), "mcmc.list")
})

test_that("the sampler's own diagnostics refuse a missing stanfit by name", {
  ds <- nosf_case()
  expect_error(nuts_params(ds), "these draws have none", fixed = TRUE)
  expect_error(log_posterior(ds), "log_posterior() reads the sampler's",
               fixed = TRUE)
})
