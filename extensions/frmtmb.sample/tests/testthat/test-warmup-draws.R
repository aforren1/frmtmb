# The warmup frm_sample() keeps (rstan's and brms's save_warmup = TRUE)
# and brms's inc_warmup on the draws accessors. Every method that does
# not ask for the warmup reads the post-warmup draws alone, which is
# asserted by comparing each one with the same object stripped of its
# warmup.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

wu_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(31)
      d <- data.frame(x = stats::rnorm(50), g = factor(rep(1:5, 10)))
      d$y <- 1 + 0.5 * d$x + stats::rnorm(5)[d$g] + stats::rnorm(50)
      run <- function(...) {
        suppressWarnings(suppressMessages(frm_sample(
          bf(y ~ x + (1 | g)), data = d, family = gaussian(),
          prior = set_prior("normal(0, 1)", class = "b"), chains = 2,
          iter = 300, warmup = 200, seed = 5, refresh = 0, ...)))
      }
      ds <- run(sample_prior = "yes")
      nw <- run(save_warmup = FALSE)
      bare <- ds
      bare$warmup_recipe <- NULL
      bare$warmup_prior <- NULL
      cache <<- list(ds = ds, nw = nw, bare = bare)
    }
    cache
  }
})

test_that("inc_warmup = TRUE puts each chain's warmup before its draws", {
  cs <- wu_case()
  ds <- cs$ds
  a <- as_draws_array(ds, inc_warmup = TRUE)
  # 200 warmup and 100 kept iterations per chain, as brms's fit has
  expect_identical(dim(a), c(300L, 2L, length(variables(ds))))
  expect_identical(nrow(as_draws_df(ds, inc_warmup = TRUE)), 600L)
  expect_identical(posterior::ndraws(as_draws_list(ds, inc_warmup = TRUE)),
                   600L)
  expect_identical(posterior::ndraws(as_draws_rvars(ds, inc_warmup = TRUE)),
                   600L)
  expect_identical(posterior::ndraws(as_draws(ds, inc_warmup = TRUE)), 600L)
  # the last 100 iterations of each chain are the draws themselves
  post <- as_draws_array(ds)
  expect_identical(unclass(a)[201:300, , ], unclass(post)[, , ],
                   ignore_attr = TRUE)
  # the first are the sampler's own warmup, moved to the draws' names:
  # b_x is a sampled coordinate the sampler and the draws share
  sf <- rstan::extract(ds$stanfit, permuted = FALSE, inc_warmup = TRUE)
  bx <- match("b_x", variables(ds))
  expect_identical(unname(unclass(a)[1:200, , bx]),
                   unname(sf[1:200, , bx]))
  # and the warmup has prior draws too, as brms's generated quantities
  expect_true(all(is.finite(unclass(a)[1:200, , "prior_b"])))
  # brms's as.data.frame() and as.matrix() pass it on
  expect_identical(nrow(as.data.frame(ds, inc_warmup = TRUE)), 600L)
  expect_identical(nrow(as.matrix(ds, inc_warmup = TRUE)), 600L)
})

test_that("as.mcmc(inc_warmup = TRUE) is brms's, refusals included", {
  cs <- wu_case()
  mc <- suppressWarnings(as.mcmc(cs$ds, inc_warmup = TRUE))
  expect_length(mc, 2L)
  expect_identical(dim(mc[[1]]), c(300L, length(variables(cs$ds))))
  expect_error(as.mcmc(cs$ds, inc_warmup = TRUE, combine_chains = TRUE),
               "Cannot include warmup samples when 'combine_chains' is TRUE.",
               fixed = TRUE)
})

test_that("without saved warmup, inc_warmup gives the draws, as brms", {
  cs <- wu_case()
  expect_null(frmtmb.sample:::draws_warmup(cs$nw))
  expect_identical(dim(as_draws_array(cs$nw, inc_warmup = TRUE)),
                   dim(as_draws_array(cs$nw)))
  expect_identical(dim(suppressWarnings(as.mcmc(cs$nw,
                                                inc_warmup = TRUE))[[1]]),
                   c(100L, length(variables(cs$nw))))
  # nsamples() counts what rstan saved
  expect_identical(suppressWarnings(nsamples(cs$nw, incl_warmup = TRUE)),
                   200L)
  expect_identical(suppressWarnings(nsamples(cs$ds, incl_warmup = TRUE)),
                   600L)
  # save_warmup = FALSE changes nothing that is kept
  expect_identical(cs$nw$draws, cs$ds$draws[, colnames(cs$nw$draws)])
})

test_that("everything else reads the post-warmup draws alone", {
  cs <- wu_case()
  ds <- cs$ds
  bare <- cs$bare
  expect_identical(ndraws(ds), 200L)
  expect_identical(niterations(ds), 100L)
  expect_identical(variables(ds), variables(bare))
  expect_identical(summary(ds), summary(bare))
  expect_identical(rhat(ds), rhat(bare))
  expect_identical(neff_ratio(ds), neff_ratio(bare))
  expect_identical(fixef(ds), fixef(bare))
  expect_identical(posterior_summary(ds), posterior_summary(bare))
  expect_identical(as_draws_df(ds), as_draws_df(bare))
  expect_identical(as.matrix(ds), as.matrix(bare))
  # bayesplot's stanfit methods take their own inc_warmup, off by default
  expect_identical(nrow(nuts_params(ds)), 2L * 100L * 6L)
  expect_identical(nrow(log_posterior(ds)), 2L * 100L)
  # and the diagnostics over the warmup are brms's, through as_draws
  rw <- rhat(ds, inc_warmup = TRUE)
  expect_identical(names(rw), names(rhat(ds)))
  a <- as_draws_array(ds, inc_warmup = TRUE)
  expect_identical(rw[["b_x"]], posterior::rhat(unclass(a)[, , "b_x"]))
})

test_that("the warmup is stored once, in the stanfit, as brms keeps it", {
  cs <- wu_case()
  ds <- cs$ds
  # no second copy in the draws' coordinates; the prior draws of its
  # rows are the one thing the stanfit cannot give back
  expect_null(ds$warmup)
  expect_identical(dim(ds$warmup_prior), c(400L, 4L))
  # the recipe is the one the draws were made with: run on the
  # post-warmup iterations it gives the draws back, bitwise
  a <- rstan::extract(ds$stanfit, permuted = FALSE)
  post <- frmtmb.sample:::sample_finish(a, ds$fit, ds$warmup_recipe)
  keep <- !startsWith(colnames(ds$draws), "prior_")
  expect_identical(post, ds$draws[, keep])
})

test_that("a warmup that does not fit the draws is refused, not dropped", {
  cs <- wu_case()
  x <- cs$ds
  # a method that edits the draws after frm_sample() leaves the warmup
  # with other variables than the draws
  x$draws <- x$draws[, -1L]
  expect_error(as_draws_array(x, inc_warmup = TRUE),
               "the saved warmup does not fit these draws")
  # and the warmup must split over the chains
  y <- cs$ds
  y$warmup <- frmtmb.sample:::draws_warmup(y)[-1L, ]
  expect_error(as.mcmc(y, inc_warmup = TRUE),
               "the saved warmup does not fit these draws")
  # without inc_warmup nothing reads it
  expect_identical(dim(as_draws_array(x)), c(100L, 2L, ncol(x$draws)))
})
