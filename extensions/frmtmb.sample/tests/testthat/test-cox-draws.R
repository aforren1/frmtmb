# A cox() baseline on the sampler (dev/moparity-findings.md, item 4).
# On 0.69.0 (rellib-r7) the baseline was sampled in the fit's softmax
# coordinates with no density on them: on seed 2 of
# dev/moparity-cox-sample.R two coordinates ran to -3e11 and -9e12, with
# 239 divergences and R-hat 3.07. It is now sampled in softmax
# coordinates under brms's dirichlet(1) and reported as brms's sbhaz[k];
# the same script compares the draws with brms 2.23.0's.

cox_draws_data <- function() {
  set.seed(2)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
}

test_that("a cox() baseline is drawn as brms's sbhaz weights", {
  skip_on_cran()
  skip_sampler()
  fit <- frm(bf(time | cens(cens) ~ x), family = cox(),
             data = cox_draws_data())
  ds <- allow_warnings(suppressMessages(frm_sample(
    fit, chains = 2, iter = 600, seed = 1, refresh = 0)),
    "Effective Samples Size")
  nm <- paste0("sbhaz[", 1:5, "]")
  expect_true(all(nm %in% variables(ds)))
  expect_false(any(grepl("^sbhaz_raw", variables(ds))))
  W <- as.matrix(ds, variable = nm)
  expect_lt(max(abs(rowSums(W) - 1)) / max(W), 1e-12)
  # every weight is interior under dirichlet(1), as in brms, where the
  # fit puts some at 0
  expect_gt(min(colMeans(W)), 0.01 * max(colMeans(W)))
  expect_lt(min(cox_baseline(fit)), 1e-6 * max(cox_baseline(fit)))
  sp <- rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)
  expect_identical(sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0)),
                   0)
  # the draws go back to the model in the fit's chart: the baseline of
  # a draw read through cox_baseline() is the draw's own weights
  im <- frmtmb.sample:::draws_internal_matrix(ds)
  z <- im[1, grep("^sbhaz", colnames(im))[1:4]]
  expect_equal(unname(frmtmb:::mo_simplex(unname(z))), unname(W[1, ]))
})

test_that("a sbhaz prior replaces dirichlet(1) on the baseline", {
  fit <- frm(bf(time | cens(cens) ~ x), family = cox(),
             data = cox_draws_data())
  pr <- set_prior("dirichlet(c(5, 1, 1, 1, 1))", class = "sbhaz")
  ent <- frmtmb::resolve_prior_input(fit, pr)$entries
  expect_identical(ent[[1L]]$comp, "sbhaz_raw")
  se <- frmtmb.sample:::mo_sample_entries(ent)
  expect_identical(se[[1L]]$dist$chart, "softmax")
  expect_null(frmtmb.sample:::mo_simplex_nlp(fit, se))
})

# The default prior note on a baseline steps aside for a sbhaz row
# (punch round 1, m4). Seen to fail on the lane build before the fix,
# which took no prior to read the row from.
test_that("the sbhaz note is left out when a sbhaz row priors the baseline", {
  fit <- frm(bf(time | cens(cens) ~ x), family = cox(),
             data = cox_draws_data())
  hit <- function(n) any(grepl("sbhaz: dirichlet(1)", n, fixed = TRUE))
  expect_true(hit(frmtmb.sample:::default_prior_notes(fit)))
  pr <- set_prior("dirichlet(2)", class = "sbhaz")
  expect_false(hit(frmtmb.sample:::default_prior_notes(fit, pr)))
})
