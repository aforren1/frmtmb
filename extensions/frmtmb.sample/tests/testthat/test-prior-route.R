# What this package contributes to frmtmb::get_prior().
#
# The registration is a load-time act with no return value, so the only
# way to see that it worked is to ask get_prior() for the route it
# answers. The companion assertion lives in frmtmb's own suite: that the
# fit route is unaffected by any registration at all.

prior_route_data <- function() {
  set.seed(7)
  data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
             g = factor(rep(1:6, 10)))
}

test_that("get_prior(route = 'sample') reports this package's defaults", {
  dd <- prior_route_data()
  form <- frmtmb::bf(y ~ x + (1 | g)) + stats::gaussian()

  gp <- frmtmb::get_prior(form, data = dd, route = "sample")
  expect_identical(attr(gp, "route"), "sample")

  # brms 2.23's defaults: a Student-t on the intercept centered on the
  # response, a half-Student-t on sigma and on the variance components
  icpt <- gp$prior[gp$class == "Intercept" & gp$dpar == "" &
                     gp$coef == ""]
  expect_match(icpt, "^student_t\\(3, ")
  sigma <- gp$prior[gp$class == "Intercept" & gp$dpar == "sigma"]
  expect_match(sigma, "^student_t\\(3, 0, ")
  expect_true(all(grepl("^student_t\\(3, 0, ",
                        gp$prior[gp$class == "sd" & gp$coef == ""])))
  # population-level slopes are flat in brms and are flat here, so the
  # column is reporting rather than filling
  expect_true(all(gp$prior[gp$class == "b"] == "(flat)"))
})

test_that("the fit route stays flat with this package loaded", {
  dd <- prior_route_data()
  form <- frmtmb::bf(y ~ x + (1 | g)) + stats::gaussian()

  gf <- frmtmb::get_prior(form, data = dd)
  expect_identical(attr(gf, "route"), "fit")
  expect_true(all(gf$prior == "(flat)"))
  # the two routes describe the same design and differ only in the
  # column that is about a route. as.data.frame() drops the label along
  # with the class, so this compares the tables rather than the routes
  # they are labeled with
  gs <- frmtmb::get_prior(form, data = dd, route = "sample")
  expect_identical(as.data.frame(gf)[, -1L], as.data.frame(gs)[, -1L])
})

test_that("a correlated block gets the LKJ default on the sample route", {
  dd <- prior_route_data()
  form <- frmtmb::bf(y ~ x + (x | g)) + stats::gaussian()
  gp <- frmtmb::get_prior(form, data = dd, route = "sample")
  expect_true(all(gp$prior[gp$class == "cor"] == "lkj(1)"))
})
