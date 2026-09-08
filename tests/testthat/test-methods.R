fit_sleep <- local({
  if (requireNamespace("lme4", quietly = TRUE)) {
    data(sleepstudy, package = "lme4")
    frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
           data = sleepstudy)
  } else {
    NULL
  }
})

test_that("accessor methods are consistent", {
  skip_if(is.null(fit_sleep))
  fit <- fit_sleep

  expect_identical(stats::nobs(fit), 180L)
  expect_s3_class(logLik(fit), "logLik")
  expect_identical(attr(logLik(fit), "df"), length(fit$opt$par))

  fe <- fixef(fit)
  expect_named(fe, c("mu", "sigma"))
  expect_named(fe$mu, c("(Intercept)", "Days"))

  V <- vcov(fit)
  expect_identical(dim(V), c(3L, 3L))
  expect_true(isSymmetric(V, tol = 1e-8))
  expect_true(all(diag(V) > 0))

  re <- ranef(fit)
  expect_length(re, 1)
  expect_identical(dim(re[[1]]), c(18L, 2L))

  vc <- VarCorr(fit)
  expect_length(vc, 1)
  expect_identical(dim(vc[[1]]), c(2L, 2L))

  expect_identical(family(fit)$family, "gaussian")
  expect_s3_class(formula(fit), "formula")
})

test_that("predict/fitted/residuals invariants hold", {
  skip_if(is.null(fit_sleep))
  fit <- fit_sleep

  mu <- predict(fit, type = "response")
  expect_length(mu, stats::nobs(fit))
  expect_identical(mu, fitted(fit))
  expect_identical(predict(fit, type = "link"), mu)  # identity link

  sig <- predict(fit, type = "response", dpar = "sigma")
  expect_true(all(sig > 0))
  expect_lt(stats::sd(sig), 1e-10)  # intercept-only sigma is constant

  r <- residuals(fit)
  expect_equal(r, fit$frame$y$Reaction - fitted(fit))
  rp <- residuals(fit, type = "pearson")
  expect_equal(rp, r / sig, tolerance = 1e-10)
})

test_that("print and summary run without error", {
  skip_if(is.null(fit_sleep))
  expect_output(print(fit_sleep), "frmtmb fit")
  s <- summary(fit_sleep)
  expect_s3_class(s, "summary.frmtmb_fit")
  expect_output(print(s), "Coefficients")
  expect_named(s$coefficients, c("mu", "sigma"))
})

test_that("start values are validated", {
  skip_if(is.null(fit_sleep))
  data(sleepstudy, package = "lme4")
  expect_error(frm(bf(Reaction ~ Days) + gaussian(), sleepstudy,
                      start = list(bogus = 1)),
               "Unknown start component")
  expect_error(frm(bf(Reaction ~ Days) + gaussian(), sleepstudy,
                      start = list(beta = 1)),
               "length")
})

# --- every ranef() block is addressable -------------------------------
# 0.52.0 re-keyed ranef() by the GROUPING FACTOR, which is brms's and
# lme4's key and stays. A nonlinear model with (1 | id) on three
# parameters then has three entries all called "id", and [["id"]]
# reached the first one silently. The key now also takes the block
# label, and the bare factor name is refused rather than answered wrong.

ranef_multi_fit <- function() {
  set.seed(9)
  w <- seq(1, 6, by = 0.5)
  d <- do.call(rbind, lapply(1:10, function(i) {
    data.frame(id = i, group = i %% 2, w = w, logw = log(w),
               I = rexp(length(w), rate = 1 / exp(1.2 - 1.5 * log(w))))
  }))
  d$id <- factor(d$id)
  frm(bf(I ~ apo - chi * logw, apo ~ 1 + (1 | id),
         chi ~ 1 + group + (1 | id), nl = TRUE),
      family = exponential(link = "log"), data = d)
}

test_that("blocks sharing a grouping factor are addressed by their term label", {
  fit <- ranef_multi_fit()
  re <- ranef(fit)
  expect_named(re, c("id", "id"))
  # the label reaches its own block, and they are different blocks
  a <- re[["apo: 1 | id"]]
  c2 <- re[["chi: 1 | id"]]
  expect_equal(dim(a), c(10L, 1L))
  expect_false(isTRUE(all.equal(as.vector(a), as.vector(c2))))
  # which is what positional indexing gives, unchanged
  expect_identical(a, re[[1]])
  expect_identical(c2, re[[2]])
  # `$` takes the label too
  expect_identical(re$`chi: 1 | id`, c2)
})

test_that("an ambiguous grouping factor is refused, not answered with the first block", {
  fit <- ranef_multi_fit()
  re <- ranef(fit)
  expect_error(re[["id"]], "2 random-effect blocks on grouping factor 'id'")
  expect_error(re[["id"]], "apo: 1 \\| id")
  expect_error(re$id, "Address a block by its term label")
})

test_that("one block per factor keeps the 0.52.0 key, and gains the label", {
  set.seed(1)
  dd <- data.frame(x = stats::rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- stats::rnorm(100, 1 + 0.5 * dd$x +
                         stats::rnorm(10, 0, 0.8)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  re <- ranef(fit)
  expect_equal(dim(re$g), c(10L, 1L))
  expect_identical(re[["g"]], re[[1]])
  expect_identical(re[["1 | g"]], re[[1]])
  # base `$` partial matching is not lost
  expect_identical(re$g, re[[1]])
  # a name that is neither is still NULL through `$`
  expect_null(re$nosuchfactor)
})

test_that("as.data.frame() and print() still index by position", {
  fit <- ranef_multi_fit()
  re <- ranef(fit)
  df <- as.data.frame(re)
  expect_equal(nrow(df), 20L)
  expect_setequal(unique(df$grp), c("apo: 1 | id", "chi: 1 | id"))
  expect_output(print(re), "apo: 1 \\| id")
})

# The multi-block case is not a nonlinear-model corner. Three ordinary
# spellings put two blocks on one grouping factor, and the ranef()
# refusal reaches all of them; only the nonlinear one was covered above.

ranef_two_block_data <- function(seed = 21, n = 200) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, n / 20)))
  d$y <- stats::rnorm(n, 1 + 0.5 * d$x +
                        stats::rnorm(20, 0, 0.8)[d$g], exp(0.2))
  d
}

test_that("an uncorrelated slope puts two blocks on one factor", {
  d <- ranef_two_block_data()
  fit <- suppressWarnings(frm(bf(y ~ x + (1 + x || g)) + gaussian(), data = d))
  re <- ranef(fit)
  expect_named(re, c("g", "g"))
  expect_identical(unname(vapply(re, function(m) attr(m, "term"), "")),
                   c("1 | g", "0 + x | g"))
  expect_error(re$g, "2 random-effect blocks on grouping factor 'g'")
  expect_identical(re[["1 | g"]], re[[1]])
  expect_identical(re[["0 + x | g"]], re[[2]])

  # the two bars written out are the same two blocks, so they behave
  # identically - this is `||` desugared, not a separate rule
  fit2 <- suppressWarnings(frm(bf(y ~ x + (1 | g) + (0 + x | g)) + gaussian(),
                               data = d))
  re2 <- ranef(fit2)
  expect_identical(unname(vapply(re2, function(m) attr(m, "term"), "")),
                   c("1 | g", "0 + x | g"))
  expect_error(re2[["g"]], "2 random-effect blocks")
})

test_that("a random effect on a second dpar puts two blocks on one factor", {
  d <- ranef_two_block_data()
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | g), sigma ~ (1 | g)) +
                                gaussian(), data = d))
  re <- ranef(fit)
  expect_named(re, c("g", "g"))
  # the dpar rides in the label, which is what tells the two apart
  expect_identical(unname(vapply(re, function(m) attr(m, "term"), "")),
                   c("1 | g", "sigma: 1 | g"))
  expect_error(re$g, "2 random-effect blocks on grouping factor 'g'")
  expect_identical(re[["sigma: 1 | g"]], re[[2]])
  expect_equal(dim(re[["sigma: 1 | g"]]), c(20L, 1L))
})

test_that("a CORRELATED slope is one block and still answers to the factor", {
  d <- ranef_two_block_data()
  fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = d))
  re <- ranef(fit)
  expect_named(re, "g")
  expect_equal(dim(re$g), c(20L, 2L))
  expect_identical(re[["1 + x | g"]], re[[1]])
})
