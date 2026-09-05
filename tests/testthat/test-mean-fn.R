# A structured family that declares no mean must not get one by accident:
# with no `post$mean_fn`, the mean is the dpar named mu by convention, and
# a family with no such dpar used to be read as "the mean is mu" anyway,
# which reported a race model's first parameter as its fitted value.

test_that("a family with neither a mean function nor a mu has no fitted value", {
  skip_on_cran()
  set.seed(4)
  d <- data.frame(x = rnorm(60))
  d$y <- rexp(60, rate = exp(-(0.2 + 0.3 * d$x)))
  fam <- frmtmb_family(
    "nomean", dpars = c("lrate", "shape"), primary_dpars = "lrate",
    links = list(lrate = "identity", shape = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dgamma(y, shape = dpars[["shape"]],
                   scale = exp(dpars[["lrate"]]) / dpars[["shape"]], log = TRUE)
    },
    init_dpars = list(lrate = function(y, aterms) log(mean(y)),
                      shape = function(y, aterms) 1)
  )
  fit <- frm(bf(y ~ x) + fam, data = d)
  expect_error(fitted(fit), "not defined for family|declares no mean")
  expect_error(predict(fit, type = "response"), "declares no mean")
  # the link-scale prediction and a dpar by name are unaffected
  expect_length(predict(fit, type = "link"), 60)
  expect_length(predict(fit, dpar = "shape"), 60)
})
