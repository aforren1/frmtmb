
# bayes_R2() and loo_R2() refuse a response named twice, as every other
# method that takes `resp` does. Before, bayes_R2() answered two
# identical rows and loo_R2() two rows with different values, from two
# runs of its bootstrap (dev/reviews/2026-10-29-release.md, finding 3).
# They keep brms's spelling of a name, which a brms pin holds them to.
test_that("bayes_R2() and loo_R2() refuse a response named twice", {
  ds <- rm_draws()
  msg <- "`resp` names the response 'y' twice"
  expect_error(bayes_R2(ds, resp = c("y", "y")), msg, fixed = TRUE)
  expect_error(rm_quiet(loo_R2(ds, resp = c("y", "y"))), msg, fixed = TRUE)
  # the guard absent: one response, both, and the order named
  expect_identical(rownames(bayes_R2(ds, resp = "y")), "R2y")
  expect_identical(rownames(bayes_R2(ds)), c("R2y", "R2y2"))
  expect_identical(rownames(bayes_R2(ds, resp = c("y2", "y"))),
                   c("R2y2", "R2y"))
  expect_identical(rownames(rm_quiet(loo_R2(ds, resp = "y2", seed = 1))),
                   "R2y2")
})
