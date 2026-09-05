# A structured family that declares no mean must not get one by accident:
# reading a missing `post$mean_fn` as "the mean is mu" reported a race
# model's drift rate as its fitted value when a dpar happened to be
# named mu, and refused only when none was.

test_that("a family without a mean function has no fitted value", {
  skip_on_cran()
  set.seed(4)
  d <- data.frame(y = rexp(60, 2), x = rnorm(60))
  fam <- custom_family(
    "nomean", dpars = c("mu", "sigma"), links = c("identity", "log"),
    lpdf = function(y, dpars, aterms) dexp(y, rate = exp(-dpars[["mu"]]) / dpars[["sigma"]], log = TRUE),
    type = "continuous"
  )
  fit <- frm(bf(y ~ x) + fam, data = d)
  # the drift-shaped parameter is called mu, which is the case that leaked
  expect_error(fitted(fit), "mean")
  expect_error(predict(fit, type = "response"), "mean")
  # the link-scale prediction is unaffected
  expect_length(predict(fit, type = "link"), 60)
})
