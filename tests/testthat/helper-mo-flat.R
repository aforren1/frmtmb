# A monotonic predictor whose middle category is never observed. With
# three categories (D = 2) the two steps enter the likelihood only as
# their sum, which is 1, so the simplex's one coordinate does not enter
# it at all: its gradient and Hessian row are exactly zero on every
# data set. test-se-check.R and test-diagnostics-ux.R read it as an
# exactly flat parameter.
#
# Before lane optima those tests read a saturated softmax coordinate of
# brms_monotonic's `ls ~ mo(income) * age` (seeds 71 and 12), whose
# Hessian row was zero only because the softmax had run its weight to
# 1e-20 or below. The simplex chart of mo_simplex() reaches a weight of
# 0 at a finite coordinate with a finite curvature, so those fits keep
# every standard error now (dev/optima-findings.md).
mo_flat_data <- function(seed = 1, n = 100) {
  set.seed(seed)
  inc <- sample(c(0L, 2L), n, TRUE)
  age <- stats::rnorm(n, mean = 40, sd = 10)
  data.frame(inc = inc, age = age,
             ls = c(30, 50, 70)[inc + 1L] + stats::rnorm(n, sd = 7))
}
