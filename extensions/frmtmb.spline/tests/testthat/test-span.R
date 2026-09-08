## Drawing a curve past a ps() knot span.
##
## A ps() basis is a partition of unity only between its frozen outer
## knots. Past them it is a partial sum that decays to zero, so a curve
## drawn there bends smoothly to whatever the rest of the nonlinear body
## gives, which is the one shape a reader will not question. frmtmb says
## so at predict(newdata = ) and, since the SPLINE-SPAN lane, at
## frm_lp_basis(newdata = ) as well.
##
## These three functions are the doors a user actually draws a curve
## through, so they are the ones the statement has to reach. Two of them
## warn once and one refuses, and the difference is deliberate: a band
## past the span is visible on the page, and a peak located past it is a
## number with a standard error beside it.

sp_span_fit <- function(seed = 4242, n = 220) {
  set.seed(seed)
  d <- data.frame(t = sort(stats::runif(n)))
  d$y <- 2 + sin(2 * pi * d$t) + stats::rnorm(n, 0, 0.25)
  fit <- frmtmb::frm(
    frmtmb::bf(y ~ lev + ps(t, k = 10, pad = 0.3), lev ~ 1, nl = TRUE),
    family = stats::gaussian(), data = d)
  list(d = d, fit = fit,
       span = fit$frame$linpreds[["y.mu"]]$ps_terms[[1]]$knot_range)
}

# every warning raised while `expr` runs, as text
sp_span_warnings <- function(expr) {
  ws <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

test_that("frm_curve() past the knot span warns once, with the span", {
  skip_on_cran()
  o <- sp_span_fit()
  g_out <- data.frame(t = seq(o$span[2] + 0.05, o$span[2] + 2,
                              length.out = 15))

  ws <- sp_span_warnings(frm_curve(o$fit, newdata = g_out,
                                   simultaneous = FALSE))
  # ONCE. The seam is read once per call but evaluates the closure more
  # than once, and before this lane it said nothing at all.
  expect_length(ws, 1L)
  expect_match(ws[[1]], "frm_curve()", fixed = TRUE)
  expect_match(ws[[1]], "outside the frozen knot span")
  # the span itself, which is the number the reader needs to act
  expect_true(grepl(format(o$span[1], digits = 4), ws[[1]], fixed = TRUE))
  expect_true(grepl(format(o$span[2], digits = 4), ws[[1]], fixed = TRUE))
  # and it counts the grid, not the coefficients
  expect_match(ws[[1]], "15 of 15", fixed = TRUE)

  # it is still a curve: the warning is a diagnostic, not a refusal
  cv <- suppressWarnings(frm_curve(o$fit, newdata = g_out,
                                   simultaneous = FALSE))
  expect_s3_class(cv, "frmtmb_curve")
  expect_equal(nrow(cv), nrow(g_out))
})

test_that("inside the span the three curve functions stay quiet", {
  skip_on_cran()
  o <- sp_span_fit()
  g_in <- data.frame(t = seq(0.05, 0.95, length.out = 15))
  expect_no_warning(frm_curve(o$fit, newdata = g_in, simultaneous = FALSE))
  expect_no_warning(frm_curve_deriv(o$fit, var = "t", newdata = g_in,
                                    simultaneous = FALSE))
  expect_no_warning(frm_curve_feature(o$fit, var = "t", type = "maximum",
                                      newdata = g_in))
  # and the feature search still finds the peak of sin(2 pi t), at 0.25
  ft <- frm_curve_feature(o$fit, var = "t", type = "maximum",
                          newdata = g_in)
  expect_equal(nrow(ft), 1L)
  expect_equal(ft$.estimate, 0.25, tolerance = 0.05)
})

test_that("frm_curve_deriv() past the span warns once and says stencil", {
  skip_on_cran()
  o <- sp_span_fit()
  g_out <- data.frame(t = seq(o$span[2] + 0.05, o$span[2] + 2,
                              length.out = 15))
  ws <- sp_span_warnings(frm_curve_deriv(o$fit, var = "t",
                                         newdata = g_out,
                                         simultaneous = FALSE))
  expect_length(ws, 1L)
  expect_match(ws[[1]], "frm_curve_deriv()", fixed = TRUE)
  # the count is over the GRID the caller passed, not over the
  # three-point stencil the design is built on: 15, not 45
  expect_match(ws[[1]], "15 of 15", fixed = TRUE)
  expect_false(grepl("45 of 45", ws[[1]], fixed = TRUE))
  expect_match(ws[[1]], "standard error of exactly zero", fixed = TRUE)
})

test_that("frm_curve_deriv() past the span works on its own defaults", {
  skip_on_cran()
  o <- sp_span_fit()
  # simultaneous = TRUE is the DEFAULT, and it is the path every other
  # past-span test here used to dodge. Past the outer knot the
  # derivative design is exactly zero, so those rows carry a standard
  # error of exactly zero, and standardizing their (also exactly zero)
  # deviation by it gave NaN and killed quantile.default() with a
  # message naming neither the span nor the function.
  g_out <- data.frame(t = seq(0.5, o$span[2] + 1, length.out = 25))
  ws <- sp_span_warnings(d <- frm_curve_deriv(o$fit, var = "t",
                                              newdata = g_out, nsim = 500,
                                              seed = 1))
  expect_length(ws, 1L)
  expect_s3_class(d, "frmtmb_curve")
  expect_equal(nrow(d), nrow(g_out))

  # the band exists, is finite, and is computed over the rows that carry
  # uncertainty; the zero-se rows get a zero-width band, which is what a
  # deterministic point deserves
  expect_true(all(is.finite(d$.crit_sim)))
  expect_gt(d$.crit_sim[1], 0)
  zero <- d$.se == 0
  expect_gt(sum(zero), 0)
  expect_equal(d$.lower_sim[zero], d$.estimate[zero])
  expect_equal(d$.upper_sim[zero], d$.estimate[zero])
  # and it is still a simultaneous band: wider than the pointwise one
  expect_gt(d$.crit_sim[1], d$.crit[1])
})

test_that("a grid with no uncertainty anywhere refuses by name", {
  skip_on_cran()
  o <- sp_span_fit()
  # every row past the outer knot: nothing left to take a maximum over,
  # so there is no simultaneous band rather than a NaN one
  far <- data.frame(t = seq(o$span[2] + 5, o$span[2] + 9, length.out = 12))
  expect_error(
    suppressWarnings(frm_curve_deriv(o$fit, var = "t", newdata = far,
                                     nsim = 200)),
    "standard error of exactly zero", fixed = TRUE)
  # and the pointwise answer is still available, as the refusal says
  expect_no_error(suppressWarnings(
    frm_curve_deriv(o$fit, var = "t", newdata = far,
                    simultaneous = FALSE)))
})

test_that("a grid ending exactly on a knot is inside the span", {
  skip_on_cran()
  o <- sp_span_fit()
  # The span check used to run on the widened difference stencil, which
  # reaches a millionth of the range past both ends of the grid. A grid
  # laid exactly on the knot span was therefore warned about by
  # frm_curve_deriv() and REFUSED by frm_curve_feature(), and narrowing
  # the grid to the span, which is what the refusal tells the user to
  # do, reproduced the refusal.
  g <- data.frame(t = seq(o$span[1], o$span[2], length.out = 25))
  expect_no_warning(frm_curve(o$fit, newdata = g, simultaneous = FALSE))
  expect_no_warning(frm_curve_deriv(o$fit, var = "t", newdata = g,
                                    simultaneous = FALSE))
  expect_no_error(frm_curve_feature(o$fit, var = "t", type = "maximum",
                                    newdata = g))
  expect_no_error(frm_curve_feature(o$fit, var = "t", type = "crossing",
                                    at = 2, newdata = g))
  # on the defaults too, which is where the stencil is widest
  expect_no_warning(frm_curve_deriv(o$fit, var = "t", newdata = g,
                                    order = 2, nsim = 200, seed = 1))

  # a hair outside, and it speaks again: the check is on the grid, not
  # loosened into a tolerance that would hide a real excursion
  w <- o$span[2] - o$span[1]
  out <- data.frame(t = seq(o$span[1] - 0.01 * w, o$span[2], length.out = 25))
  expect_warning(frm_curve_deriv(o$fit, var = "t", newdata = out,
                                 simultaneous = FALSE),
                 "outside the frozen knot span")
  expect_error(frm_curve_feature(o$fit, var = "t", type = "maximum",
                                 newdata = out),
               "the search bracket leaves", fixed = TRUE)
})

test_that("frm_curve_feature() refuses a bracket that leaves the span", {
  skip_on_cran()
  o <- sp_span_fit()
  g_out <- data.frame(t = seq(o$span[2] + 0.05, o$span[2] + 2,
                              length.out = 15))
  # by name, and before any root is refined
  expect_error(frm_curve_feature(o$fit, var = "t", type = "maximum",
                                 newdata = g_out),
               "frm_curve_feature(): the search bracket leaves",
               fixed = TRUE)
  expect_error(frm_curve_feature(o$fit, var = "t", type = "crossing",
                                 at = 2, newdata = g_out),
               "outside the frozen knot span")
  # the refusal carries the span, so the remedy is in the message
  e <- tryCatch(frm_curve_feature(o$fit, var = "t", type = "maximum",
                                  newdata = g_out),
                error = function(e) conditionMessage(e))
  expect_true(grepl(format(o$span[2], digits = 4), e, fixed = TRUE))
  expect_match(e, "Narrow newdata to the span", fixed = TRUE)

  # a grid that STRADDLES the span is refused too: half a bracket past
  # the knots is still a root of the decaying partial sum
  g_half <- data.frame(t = seq(0.5, o$span[2] + 1, length.out = 15))
  expect_error(frm_curve_feature(o$fit, var = "t", type = "maximum",
                                 newdata = g_half),
               "the search bracket leaves", fixed = TRUE)
})

test_that("the refusal replaces a warning per Newton step, not adds to it", {
  skip_on_cran()
  o <- sp_span_fit()
  g_out <- data.frame(t = seq(o$span[2] + 0.05, o$span[2] + 2,
                              length.out = 15))
  # Before this lane the feature search reached the curve through
  # predict(newdata = ), which IS armed, so one call raised the span
  # warning eleven times: once for the grid scan and once per Newton
  # iteration per root. It now refuses, and raises none.
  ws <- character(0)
  # named, not bare: a bare expect_error() would pass on an error raised
  # before the span check ever ran, which is the failure this test is
  # meant to catch
  expect_error(withCallingHandlers(
    frm_curve_feature(o$fit, var = "t", type = "maximum",
                      newdata = g_out),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    "the search bracket leaves", fixed = TRUE)
  expect_length(ws, 0L)
})

test_that("a model with no ps() term is untouched by any of this", {
  skip_on_cran()
  set.seed(11)
  d <- data.frame(x = sort(stats::runif(220)))
  d$y <- 2 * sin(pi * d$x) + stats::rnorm(220, 0, 0.35)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                     family = stats::gaussian(), data = d)
  # a grid far outside the DATA range, which for an ordinary smooth is
  # extrapolation and not a partition-of-unity cliff: nothing here has
  # anything to say about it
  g <- data.frame(x = seq(5, 9, length.out = 15))
  expect_no_warning(cv <- frm_curve(fit, newdata = g, simultaneous = FALSE))
  expect_no_warning(frm_curve_deriv(fit, var = "x", newdata = g,
                                    simultaneous = FALSE))
  expect_no_error(frm_curve_feature(fit, var = "x", type = "maximum",
                                    newdata = g))
  # and the covariance check is unchanged, which is what says the
  # collector did not disturb the seam
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  expect_equal(attr(cv, "check")$n_predict, 1L)
})

## A SECOND ps() term makes the stencil guard reachable.
##
## frm_curve_feature() checks the span twice, and the two checks are not
## the same question. The scan and the five-point stencil are built from
## row 1 replicated, so they hold every column but `var` pinned; the
## re-ask runs one predict() on the WHOLE grid. A second ps() term can
## therefore leave its span in a row the scan never evaluates: the first
## refusal's guard is false, and only the stencil guard is left between
## the user and a root reported with a standard error. This case was
## deleted once as unreachable, which made the call silent where it had
## refused, so it is pinned here.
sp_span_fit2 <- function(seed = 4242, n = 220) {
  set.seed(seed)
  d <- data.frame(t = sort(stats::runif(n)), z = stats::runif(n, 0, 1))
  d$y <- 2 + sin(2 * pi * d$t) + 0.5 * d$z + stats::rnorm(n, 0, 0.25)
  fit <- frmtmb::frm(
    frmtmb::bf(y ~ lev + ps(t, k = 10, pad = 0.3) + ps(z, k = 8, pad = 0.02),
               lev ~ 1, nl = TRUE),
    family = stats::gaussian(), data = d)
  pt <- fit$frame$linpreds[["y.mu"]]$ps_terms
  list(fit = fit, t_span = pt[[1]]$knot_range, z_span = pt[[2]]$knot_range)
}

test_that("a second ps() term's span reaches the stencil refusal", {
  skip_on_cran()
  o <- sp_span_fit2()

  # The grid ends just inside t's span, so the crossing scan (which
  # reads the grid itself) stays clean, while a root within e2 of that
  # end pushes the five-point stencil past it. Default eps throughout.
  a <- 0.05
  b <- o$t_span[2] - 5e-6 * (o$t_span[2] - a)
  g <- data.frame(t = seq(a, b, length.out = 15))
  g$z <- 0.5
  expect_true(all(g$t >= o$t_span[1] & g$t <= o$t_span[2]))

  eta_at <- function(tv) {
    dd <- g[rep(1L, length(tv)), , drop = FALSE]
    dd$t <- tv
    as.numeric(stats::predict(o$fit, newdata = dd, type = "link"))
  }
  e2 <- 1e-4 * diff(range(g$t))
  at <- eta_at(b - 0.2 * e2)

  # Control. The stencil DOES leave t's span here, so the guard is true,
  # but the grid itself is clean and the re-ask correctly declines to
  # refuse. This is the half of the block that must stay quiet.
  expect_no_error(ok <- frm_curve_feature(o$fit, var = "t",
                                          type = "crossing", at = at,
                                          newdata = g))
  expect_equal(nrow(ok), 1L)

  # One row leaves z's span, and it is not row 1, so nothing the scan
  # evaluates ever sees it. Only the stencil guard can refuse here, so
  # this errors if and only if that block is present.
  g_out <- g
  g_out$z[10] <- o$z_span[2] + 5
  expect_error(frm_curve_feature(o$fit, var = "t", type = "crossing",
                                 at = at, newdata = g_out),
               "the search bracket leaves", fixed = TRUE)
})
