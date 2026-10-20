# drop1() and update() on a categorical(refcat = NA) fit with a
# transformed term. The degrees of freedom are lane catmulti's (the
# reference twin's count, taken when the fit is made, loglik_df()); the
# refits are lane drop1raw's (the data the fit holds, with its bases and
# NA rows). Neither lane could run this alone: on catmulti's build the
# refit read the model frame, which has `log(x2)` and no `x2`; on
# drop1raw's build refcat = NA counted the common shift as degrees of
# freedom. Added at the 0.72.0 consolidation (dev/round-20261020.md).

rr_cat_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), x2 = exp(stats::rnorm(n)),
                   x3 = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x - 0.4 * log(dd$x2),
             -0.4 - 0.5 * dd$x + 0.3 * log(dd$x2))
  P <- exp(E) / rowSums(exp(E))
  dd$y <- factor(apply(P, 1, function(p) {
    sample(c("a", "b", "c"), 1, prob = p)
  }), levels = c("a", "b", "c"))
  dd
}

rr_rel <- function(a, b) max(abs(a - b) / abs(b))

test_that("drop1() on categorical(refcat = NA) refits log(x2) from x2", {
  se_na <- "Standard errors are not available"
  dd <- rr_cat_data(31)
  allow_warnings(
    fn <- frm(bf(y ~ x + log(x2)), family = categorical(refcat = NA),
              data = dd),
    se_na, require = "common shift of every category's predictor")
  fa <- frm(bf(y ~ x + log(x2)), family = categorical(), data = dd)
  tn <- allow_warnings(drop1(fn, test = "Chisq"), se_na)
  ta <- drop1(fa, test = "Chisq")
  expect_identical(rownames(tn), c("<none>", "x", "log(x2)"))
  # one parameter per non-reference category, as with a reference
  expect_identical(tn$Df, ta$Df)
  expect_identical(ta$Df[-1L], c(2, 2))
  expect_lt(rr_rel(tn$AIC, ta$AIC), 1e-6)
  expect_lt(rr_rel(tn$LRT[-1L], ta$LRT[-1L]), 1e-6)
})

test_that("a refit on the fit's rows counts the twin on those rows", {
  # rows na.action dropped, and a spline basis computed on every row:
  # the refit and its twin take the fit's bases and rows, so the
  # drop1() row equals update()'s submodel
  se_na <- "Standard errors are not available"
  dd <- rr_cat_data(32)
  dd$x3[c(3, 50, 77)] <- NA
  fm <- y ~ splines::ns(x2, 3) + x3
  allow_warnings(
    fn <- frm(bf(fm), family = categorical(refcat = NA), data = dd),
    se_na)
  fa <- frm(bf(fm), family = categorical(), data = dd)
  expect_equal(nobs(fn), 297)
  tn <- allow_warnings(drop1(fn, scope = "x3"), se_na)
  ta <- drop1(fa, scope = "x3")
  expect_identical(tn$Df, ta$Df)
  expect_lt(rr_rel(tn$AIC, ta$AIC), 1e-6)
  un <- allow_warnings(update(fn, . ~ . - x3), se_na)
  ua <- update(fa, . ~ . - x3)
  expect_equal(nobs(un), 297)
  expect_identical(attr(logLik(un), "df"), attr(logLik(ua), "df"))
  expect_lt(rr_rel(AIC(un), tn$AIC[2L]), 1e-6)
  expect_lt(rr_rel(AIC(ua), ta$AIC[2L]), 1e-6)
})

test_that("drop1() stops with a refit's reason when its df is NA", {
  # no model reaches this branch today (a refit has its fit's family and
  # formulas); lane drop1raw's review asked for it to be seen to fire
  dd <- rr_cat_data(33)
  fa <- frm(bf(y ~ x + log(x2)), family = categorical(), data = dd)
  real <- loglik_df
  local_mocked_bindings(loglik_df = function(object) {
    if (identical(object, fa)) return(real(object))
    structure(NA_integer_, df_why = "a reason from the refit")
  })
  expect_error(drop1(fa),
               paste("the refit without x has no degrees of freedom to",
                     "compare: a reason from the refit"), fixed = TRUE)
})
