# fitted()'s finite-difference Est.Error on a category distribution:
# the chain rule through the linear predictor for a smooth or gp()
# block, which costs one pair of evaluations for the block where one
# pair per coefficient was the whole cost, and the kriging variance at an
# unseen gp() position, which the route left out. The cost cells and
# their counts are dev/gpby-fdcost.R's; the counts on 0.67.0 were 330
# for the 160-coefficient gp() at three new rows.

fdc_fit <- local({
  memo <- NULL
  function() {
    if (!is.null(memo)) return(memo)
    set.seed(41)
    n <- 120
    d <- data.frame(x = stats::runif(n, -2, 2))
    lat <- 3 * sin(1.5 * d$x) + stats::rlogis(n)
    d$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
                  ordered = TRUE)
    memo <<- list(d = d, fit = suppressWarnings(
      frm(bf(y ~ gp(x)), family = cumulative(), data = d)))
    memo
  }
})

# fitted()'s calls of the fitted value, the estimate's own included
fdc_count <- function(expr) {
  ns <- asNamespace("frmtmb")
  n <- 0L
  suppressMessages(trace("fitted_point", tracer = function() n <<- n + 1L,
                         where = ns, print = FALSE))
  on.exit(suppressMessages(untrace("fitted_point", where = ns)))
  force(expr)
  n
}

# The route without the chain rule and without the extra variance: one
# pair per outer parameter and per group effect. At fitted()'s own step
# (1e-5) its truncation error, which falls as the step squared, reached
# 1.2e-5 relative on a gp(x) + (1 | g) fixture, where the chain route
# moves by 3e-9 across steps from 1e-3 to 1e-6
# (dev/reviews/2026-10-05-gpby.md, m10). At a step of 1e-7 its
# truncation is about 1e-9 and its rounding about eps / 1e-7 = 2e-9, so
# it is the accurate route, and a bound of 1e-7 sits fifty times above
# both while staying far below any defect the tests pin
fdc_plain <- function(fit, nd, re_formula, b_idx = smooth_b_idx(fit)) {
  f <- function(x) fitted_point(x, nd, re_formula)
  fit_fd_se(fit, f, eps = 1e-7, b_idx = b_idx)
}
fdc_tol <- 1e-7

test_that("a gp() block off its positions costs one pair, not one each", {
  o <- fdc_fit()
  nd <- data.frame(x = c(-1.5, 0.05, 1.5))
  p <- length(outer_par_names(o$fit))
  n_b <- length(o$fit$estimates$b)
  k <- fdc_count(fv <- fitted(o$fit, newdata = nd, re_formula = NA))
  # the estimate and its reference value, a pair per outer parameter,
  # and ONE pair for the gp() block's eta; on 0.67.0 a pair per
  # coefficient besides
  expect_identical(k, 2L + 2L * p + 2L)
  expect_gt(n_b, 50L)
  expect_true(all(is.finite(fv[, "Est.Error", ])))
})

test_that("the chain rule gives the per-coefficient derivative", {
  o <- fdc_fit()
  # at observed positions there is no kriging variance, so the chain
  # route and the per-coefficient route difference the same function
  nd <- o$d[c(1, 5, 9), "x", drop = FALSE]
  fv <- fitted(o$fit, newdata = nd, re_formula = NA)
  ref <- fdc_plain(o$fit, nd, NA)
  expect_lt(max(abs(fv[, "Est.Error", ] / ref - 1)), fdc_tol)
  # and a Hilbert-space gp() in sample, every row loading every column
  fh <- suppressWarnings(frm(bf(y ~ gp(x, k = 12)), family = cumulative(),
                             data = o$d))
  k <- fdc_count(fh_fv <- fitted(fh, re_formula = NA))
  expect_identical(k, 2L + 2L * length(outer_par_names(fh)) + 2L)
  ref <- fdc_plain(fh, NULL, NA)
  expect_lt(max(abs(fh_fv[, "Est.Error", ] / ref - 1)), fdc_tol)
  # and the reviewer's fixture, a gp() beside a group term at a seen
  # level, where the step-1e-5 reference was 1.2e-5 away
  set.seed(11)
  n <- 120
  d <- data.frame(x = round(stats::runif(n, 0, 5), 1), z = stats::rnorm(n),
                  g = factor(sample(letters[1:8], n, TRUE)))
  lat <- sin(d$x) + 0.5 * d$z + stats::rnorm(8)[as.integer(d$g)] * 0.5 +
    stats::rlogis(n)
  d$y <- factor(cut(lat, stats::quantile(lat, c(0, 0.3, 0.6, 1)),
                    include.lowest = TRUE, labels = FALSE), ordered = TRUE)
  fg <- suppressWarnings(frm(bf(y ~ gp(x) + (1 | g)), data = d,
                             family = cumulative()))
  ndg <- data.frame(x = c(5.6, 6.3, 2.55),
                    g = factor(c("a", "b", "c"), levels = levels(d$g)))
  fvg <- fitted(fg, newdata = ndg)
  bi <- sort(unique(c(smooth_b_idx(fg), re_governed_b(fg))))
  # the plain route has no kriging term, so compare at the coefficient
  # part: Est.Error^2 less (dP/deta)^2 * extra_var is that part
  refg <- fdc_plain(fg, ndg, NULL, bi)
  ev <- frm_lp_basis(fg, newdata = ndg)$extra_var
  eta <- as.numeric(frm_linpred(fg, newdata = ndg, type = "link"))
  fx <- fixef(fg)
  tau <- fx[grepl("^Intercept", rownames(fx)), "Estimate"]
  gq <- cbind(-stats::dlogis(tau[1] - eta),
              stats::dlogis(tau[1] - eta) - stats::dlogis(tau[2] - eta),
              stats::dlogis(tau[2] - eta))
  coef_part <- sqrt(pmax(fvg[, "Est.Error", ]^2 - gq^2 * ev, 0))
  expect_lt(max(abs(coef_part / refg - 1)), 1e-6)
})

test_that("Est.Error carries the kriging variance past the positions", {
  o <- fdc_fit()
  # past the data the gp()'s conditional variance is most of the field's
  # own, and a category probability moves with it through eta
  nd <- data.frame(x = c(2.6, 3.2))
  fv <- fitted(o$fit, newdata = nd, re_formula = NA)
  plain <- fdc_plain(o$fit, nd, NA)
  ev <- frm_lp_basis(o$fit, newdata = nd, re_formula = NA)$extra_var
  expect_true(all(ev > 0))
  # cumulative logit: P(Y = k) = F(tau_k - eta) - F(tau_{k-1} - eta)
  eta <- as.numeric(frm_linpred(o$fit, newdata = nd, re_formula = NA,
                                type = "link"))
  fx <- fixef(o$fit)
  tau <- fx[grepl("^Intercept", rownames(fx)), "Estimate"]
  dF <- function(z) stats::dlogis(z)
  g <- cbind(-dF(tau[1] - eta),
             dF(tau[1] - eta) - dF(tau[2] - eta),
             dF(tau[2] - eta))
  want <- sqrt(plain^2 + g^2 * ev)
  expect_lt(max(abs(fv[, "Est.Error", ] / want - 1)), 1e-5)
  # the share is not small here, so leaving it out is visible
  expect_true(all(fv[, "Est.Error", ] > plain))
  expect_true(all(fv[2L, "Est.Error", ] > 1.05 * plain[2L, ]))
})
