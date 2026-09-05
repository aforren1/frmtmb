# The variable collector for a nonlinear body. `all.vars()` drops the
# ARGUMENTS of a call that sits in function position, so a variable
# named only there was never asked of `data` and the body then failed on
# it. See dev/spline-seam-proposal.md, Part 3.

test_that("nl_body_vars() sees through a call in function position", {
  bv <- frmtmb:::nl_body_vars
  expect_identical(all.vars(quote(f(x)(y))), "y")
  expect_setequal(bv(quote(f(x)(y))), c("x", "y"))
  expect_setequal(bv(quote(a * curry(tv)(zv))), c("a", "tv", "zv"))
  expect_setequal(bv(quote(g(h(p))(q)(r))), c("p", "q", "r"))
})

test_that("nl_body_vars() agrees with all.vars() on ordinary bodies", {
  bv <- frmtmb:::nl_body_vars
  bodies <- list(
    quote(a * x + b),
    quote(a * exp(-b * t)),
    quote(ult * (1 - exp(-exp(lrc) * x))),
    quote(int + exp(amp) * sin(2 * pi * (age + shift))),
    quote(ls + th * log(abs(mu))),
    quote(b0 + b1 * I(x^2)),
    quote(f(x, y = z)),
    quote(sum(w * v) / n)
  )
  for (b in bodies) {
    expect_identical(bv(b), all.vars(b), label = deparse1(b))
  }
})

test_that("nl_body_vars() keeps all.vars()'s conventions on $, @ and ::", {
  bv <- frmtmb:::nl_body_vars
  # `all.vars()` collects the field name of `$` and `@` and does not
  # collect a package name; both conventions are kept, so no body that
  # fitted before collects a different set now
  expect_identical(bv(quote(a * d$col)), all.vars(quote(a * d$col)))
  expect_identical(bv(quote(a * stats::rnorm(n))),
                   all.vars(quote(a * stats::rnorm(n))))
  expect_identical(bv(quote(a * obj@slot)), all.vars(quote(a * obj@slot)))
  # a function literal has a pairlist for formals, not a call
  expect_identical(bv(quote(sapply(x, function(u) u * k))),
                   all.vars(quote(sapply(x, function(u) u * k))))
})

test_that("nl_body_vars() survives an empty argument", {
  bv <- frmtmb:::nl_body_vars
  # `m[, 1]` holds the missing-argument symbol, which can be extracted
  # and not passed on: the error fires at the callee, so guarding the
  # extraction alone is not enough
  expect_identical(bv(quote(b0 * m[, 1])), all.vars(quote(b0 * m[, 1])))
  expect_identical(bv(quote(x[1, ])), all.vars(quote(x[1, ])))
  expect_identical(bv(quote(a * z[, , 2])), all.vars(quote(a * z[, , 2])))
  # and through a function-position call, which is the case the walker
  # exists for
  expect_setequal(bv(quote(f(m[, 1])(v))), c("m", "v"))
})

test_that("a matrix column in a nonlinear body still reaches the frame", {
  skip_on_cran()
  set.seed(4042)
  d <- data.frame(y = rnorm(50))
  m <- matrix(rnorm(2 * nrow(d)), nrow(d), 2)
  fr <- frm(bf(y ~ b0 * m[, 1], b0 ~ 1, nl = TRUE), d, gaussian(),
            dry_run = "frame", start = list(beta = 1))
  expect_true("m" %in% names(fr$linpreds[["y.mu"]]$data_list))
})

test_that("a curried call in a nonlinear body fits", {
  skip_on_cran()
  curry <- function(u) function(v) u * v
  curry2 <- function(u, v) u * v
  set.seed(4041)
  d <- data.frame(tv = runif(60), zv = runif(60))
  d$y <- 2 * d$tv * d$zv + rnorm(60, 0, 0.05)

  f_plain <- frm(bf(y ~ a * curry2(tv, zv), a ~ 1, nl = TRUE), d, gaussian())
  f_curry <- frm(bf(y ~ a * curry(tv)(zv), a ~ 1, nl = TRUE), d, gaussian())
  # identical model, so identical coefficient: before the fix the second
  # spelling died with "object 'tv' not found"
  expect_equal(fixef(f_curry)$a[[1]], fixef(f_plain)$a[[1]], tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f_curry)), as.numeric(logLik(f_plain)),
               tolerance = 1e-8)
})
