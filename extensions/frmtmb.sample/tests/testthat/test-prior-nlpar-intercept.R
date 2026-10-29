# frmtmb's spelling of a nonlinear parameter's intercept,
# set_prior(class = "Intercept", nlpar = "a"), on frm_sample(). Where
# a ~ 1 it is brms's class row class = "b", nlpar = "a", bound
# included, so the bound reaches the prior draws as well as the
# posterior. Through the first fix round it was a coefficient row whose
# bound the prior draws ignored: the prior of a went down to -1.67
# under lb = 0 (dev/reviews/2026-10-28-meme.md, S1). Where the
# parameter has more than an intercept, a bound on its intercept is a
# bound on one coefficient, which brms refuses, and so does this.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

pni_data <- function() {
  set.seed(2)
  d <- data.frame(x = stats::runif(60, 0, 3), w = stats::rnorm(60))
  d$y <- 2 * exp(-0.7 * d$x) + stats::rnorm(60, 0, 0.1)
  d
}

pni_form <- function() bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)

pni_ri <- function(fit, prior) {
  suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, prior))$ri
}

pni_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("on a ~ 1 the intercept spelling is brms's class row", {
  d <- pni_data()
  fit <- frm(pni_form(), data = d, dry_run = "objective")
  pf <- set_prior("normal(1, 1)", class = "Intercept", nlpar = "a",
                  lb = 0) +
    set_prior("normal(1, 1)", class = "Intercept", nlpar = "b", lb = 0.4)
  pb <- set_prior("normal(1, 1)", class = "b", nlpar = "a", lb = 0) +
    set_prior("normal(1, 1)", class = "b", nlpar = "b", lb = 0.4)
  rf <- pni_ri(fit, pf)
  rb <- pni_ri(fit, pb)
  # the same bounds and densities as brms's spelling (the guard-absent
  # case: class = "b", nlpar = is unchanged)
  expect_identical(rf$lower, rb$lower)
  expect_identical(rf$upper, rb$upper)
  expect_identical(lapply(rf$entries, `[[`, "dist"),
                   lapply(rb$entries, `[[`, "dist"))
  expect_identical(unname(rb$lower), c(0, 0.4))
  # the table the prior draws are read from carries the bound on the
  # class row, where brms draws prior_b_a under it
  sp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pf))
  tab <- validate_prior(sp$effective, fit)
  row <- tab[tab$class == "b" & tab$nlpar == "a" & !nzchar(tab$coef), ]
  expect_identical(as.numeric(row$lb), 0)
  expect_identical(row$prior, "normal(1, 1)")
})

test_that("the prior draws and the posterior respect the bound", {
  skip_sampler()
  d <- pni_data()
  ds <- pni_quiet(frm_sample(
    pni_form(), data = d, family = gaussian(),
    prior = set_prior("normal(1, 1)", class = "Intercept", nlpar = "a",
                      lb = 0) +
      set_prior("normal(1, 1)", class = "Intercept", nlpar = "b",
                lb = 0.4) +
      set_prior("exponential(1)", class = "sigma"),
    sample_prior = "yes", chains = 2, iter = 600, seed = 3, refresh = 0))
  pd <- prior_draws(ds)
  # brms's names for a class row's prior draws
  expect_true(all(c("b_a", "b_b") %in% names(pd)))
  expect_gte(min(pd$b_a), 0)
  expect_gte(min(pd$b_b), 0.4)
  m <- as.matrix(ds)
  expect_gte(min(m[, "b_a_Intercept"]), 0)
  expect_gte(min(m[, "b_b_Intercept"]), 0.4)
  # each draw is its normal(1, 1) truncated at the bound
  tn <- function(lo) {
    function(q) {
      (stats::pnorm(q, 1, 1) - stats::pnorm(lo, 1, 1)) /
        stats::pnorm(lo, 1, 1, lower.tail = FALSE)
    }
  }
  ks <- function(v, cdf) suppressWarnings(stats::ks.test(v, cdf))$p.value
  expect_gt(ks(pd$b_a, tn(0)), 0.01)
  expect_gt(ks(pd$b_b, tn(0.4)), 0.01)
})

test_that("with a slope, a bound on the intercept is refused", {
  d <- pni_data()
  f <- bf(y ~ a * exp(-b * x), a ~ 1 + w, b ~ 1, nl = TRUE)
  fit <- frm(f, data = d, dry_run = "objective")
  msg <- "Prior argument 'coef' may not be specified when using boundaries"
  expect_error(pni_ri(fit, set_prior("normal(1, 1)", class = "Intercept",
                                     nlpar = "a", lb = 0)),
               msg, fixed = TRUE)
  # brms's own spelling of that row is refused the same way
  expect_error(pni_ri(fit, set_prior("normal(1, 1)", class = "b",
                                     coef = "Intercept", nlpar = "a",
                                     lb = 0)),
               msg, fixed = TRUE)
  # b ~ 1 beside it is still a class row, and a's class row takes a bound
  ri <- pni_ri(fit, set_prior("normal(1, 1)", class = "Intercept",
                              nlpar = "b", lb = 0.4) +
                 set_prior("normal(0, 1)", class = "b", nlpar = "a",
                           lb = -1))
  expect_identical(unname(ri$lower[grepl("^b_", names(ri$lower))]), 0.4)
  # without a bound the intercept's own density is a coefficient row
  ri2 <- pni_ri(fit, set_prior("normal(1, 1)", class = "Intercept",
                               nlpar = "a"))
  expect_length(Filter(function(e) identical(e$dist$kind, "normal"),
                       ri2$entries), 1L)
})
