# A list of families as `family =`, one per response, as brms takes it
# for a multivariate model (brms:::validate_formula.mvbrmsformula()).
# Each entry fills its response when that response has no family, as a
# single family argument does.

family_list_data <- function() {
  set.seed(8)
  d <- data.frame(x = rnorm(60))
  d$y1 <- rnorm(60, 1 + 0.5 * d$x)
  d$y2 <- rpois(60, exp(0.3 + 0.2 * d$x))
  d
}

test_that("a list of families gives each response its own", {
  d <- family_list_data()
  fit_list <- frm(bf(y1 ~ x) + bf(y2 ~ x), data = d,
                  family = list(gaussian, poisson()))
  fit_plus <- frm(bf(y1 ~ x) + gaussian() + bf(y2 ~ x) + poisson(),
                  data = d)
  expect_identical(vapply(fit_list$spec$responses,
                          function(r) r$family$family, ""),
                   c(y1 = "gaussian", y2 = "poisson"))
  expect_identical(as.numeric(logLik(fit_list)),
                   as.numeric(logLik(fit_plus)))
  # an entry fills only a response that has no family of its own
  sp <- frm(bf(y1 ~ x) + bf(y2 ~ x) + poisson(), data = d,
            family = list(gaussian(), gaussian()), dry_run = "spec")
  expect_identical(vapply(sp$responses, function(r) r$family$family, ""),
                   c(y1 = "poisson", y2 = "poisson"))
})

test_that("the prior table of a family list is brms's", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- family_list_data()
  fl <- list(gaussian, poisson())
  pf <- get_prior(bf(y1 ~ x) + bf(y2 ~ x), data = d, family = fl)
  pb <- suppressMessages(brms::default_prior(
    brms::bf(y1 ~ x) + brms::bf(y2 ~ x), data = d, family = fl))
  key <- function(p) sort(paste(p$class, p$coef, p$resp, p$dpar))
  expect_identical(key(pf), key(as.data.frame(pb)))
})

test_that("a family list that does not fit the formula is refused", {
  d <- family_list_data()
  expect_error(frm(bf(y1 ~ x) + bf(y2 ~ x), data = d,
                   family = list(gaussian())),
               "it has to be of the same length as the number of response",
               class = "frmtmb_error")
  expect_error(frm(y1 ~ x, data = d, family = list(gaussian())),
               "A list of families is for a multivariate model",
               class = "frmtmb_error")
  expect_error(frm(bf(y1 ~ x) + bf(y2 ~ x), data = d,
                   family = list(gaussian(), 3)),
               "Cannot interpret `family`", class = "frmtmb_error")
})
