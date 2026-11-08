# Review of lane arstep: a test the lane's file lacks. A second descent that
# LOSES must leave the objective where the first descent ended, or the
# modes, fitted values and every reader built from the objective's last
# point read the second descent's end. The lane's test 5 checks the state
# only on a row where the second descent WINS, where the mutant that drops
# obj_settle() changes nothing (dev/review-arstep-out/statecmp.txt).
# Run by dev/review-arstep-mut.R <mutant> <this file>.

test_that("a second descent that loses leaves the modes at the kept optimum", {
  skip_on_cran()
  # the ncgp review's Gamma row (seed 1, n 40, curve sin(x) + 0.4 sin(10x)):
  # the working marginal's second maximum is 18.8 below its best, and the
  # descent from it ends about 20 below the first
  n <- 40L
  set.seed(74000 + 31 + n)
  x <- sort(round(stats::runif(n, 0, 10), 4))
  eta <- sin(x) + 0.4 * sin(10 * x)
  d <- data.frame(x = x, y = stats::rgamma(n, 400, 400 / exp(eta)))
  fit <- frm(y ~ gp(x), data = d, family = Gamma(link = "log"))
  # the modes at the reported optimum, solved afresh from the same point
  ref <- frm(y ~ gp(x), data = d, family = Gamma(link = "log"),
             start = list(theta = fit$estimates$theta,
                          beta = fit$estimates$beta))
  b <- fit$estimates$b
  br <- ref$estimates$b
  # whether or not the second descent ran or won here, the modes and the
  # fitted values must be those of the optimum the fit reports
  expect_lt(max(abs(b - br)), 1e-4 * max(abs(br)))
  fe <- function(f) {
    v <- fitted(f)
    as.numeric(if (is.matrix(v) || is.data.frame(v)) v[, 1] else v)
  }
  f1 <- fe(fit)
  f2 <- fe(ref)
  expect_lt(max(abs(f1 - f2)), 1e-4 * max(abs(f2)))
  # and the objective's state names that optimum
  env <- fit$obj$env
  expect_equal(unname(env$last.par[-env$random]), unname(fit$opt$par))
})
