# conditional_effects() options brms has: surface, too_far,
# select_points, spaghetti, and the Wald band of a nonlinear predictor.

ceo_data <- function(n = 150, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  f = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- 1 + d$x - d$z + 0.5 * (d$f == "b") + stats::rnorm(n, 0, 0.3)
  d$y2 <- 2 * exp(0.8 * d$z) * (1 + 0.3 * (d$f == "b")) +
    stats::rnorm(n, 0, 0.3)
  d
}

test_that("surface = TRUE draws both predictors over their range", {
  d <- ceo_data()
  fit <- frm(bf(y ~ x * z), family = gaussian(), data = d)
  ce <- conditional_effects(fit, "x:z", surface = TRUE, resolution = 9)
  df <- ce[["x:z"]]
  expect_true(attr(df, "surface"))
  expect_equal(nrow(df), 81)
  expect_true(is.numeric(df$effect2__))
  expect_equal(range(df$z), range(d$z))
  nd <- df
  expect_equal(df$estimate__, as.vector(frm_linpred(fit, newdata = nd,
                                                    type = "response")))
  # too_far removes exactly the points mgcv calls too far
  far <- mgcv::exclude.too.far(df$x, df$z, d$x, d$z, dist = 0.15)
  ct <- conditional_effects(fit, "x:z", surface = TRUE, resolution = 9,
                            too_far = 0.15)[["x:z"]]
  expect_equal(nrow(ct), sum(!far))
  expect_equal(ct$estimate__, df$estimate__[!far])
  # without surface, too_far changes nothing, as in brms
  expect_identical(
    conditional_effects(fit, "x:z", resolution = 9, too_far = 0.15),
    conditional_effects(fit, "x:z", resolution = 9))
  expect_false(attr(conditional_effects(fit, "x:z")[[1]], "surface"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(ce, ask = FALSE))
})

test_that("select_points keeps the observations near the held values", {
  d <- ceo_data()
  fit <- frm(bf(y ~ f + x + z), family = gaussian(), data = d)
  all_pts <- attr(conditional_effects(fit, "f")[[1]], "points")
  expect_equal(nrow(all_pts), nrow(d))
  expect_identical(names(all_pts), c("f", "resp__", "effect1__"))
  sel <- attr(conditional_effects(fit, "f", select_points = 0.2)[[1]],
              "points")
  unit <- function(v, at) abs((v - min(v)) / diff(range(v)) -
                                (at - min(v)) / diff(range(v)))
  keep <- unit(d$x, mean(d$x)) <= 0.2 & unit(d$z, mean(d$z)) <= 0.2
  expect_equal(nrow(sel), sum(keep))
  expect_lt(nrow(sel), nrow(d))
  expect_equal(sort(sel$resp__), sort(d$y[keep]))
  expect_error(conditional_effects(fit, "f", select_points = -1),
               "`select_points` must be one number")
})

test_that("spaghetti needs draws: boot refits have them, Wald does not", {
  d <- ceo_data()
  fit <- frm(bf(y ~ x + z), family = gaussian(), data = d)
  expect_error(conditional_effects(fit, "x", spaghetti = TRUE),
               "band = \"boot\"")
  expect_error(conditional_effects(fit, "x", spaghetti = TRUE,
                                   surface = TRUE),
               "Cannot use 'spaghetti' and 'surface' at the same time")
  ce <- conditional_effects(fit, "x", spaghetti = TRUE, band = "boot",
                            boot = 6, seed = 1, resolution = 5)
  sp <- attr(ce$x, "spaghetti")
  expect_equal(nrow(sp), 6 * 5)
  expect_identical(levels(sp$sample__), as.character(1:6))
  # each curve is one refit's: the percentiles of the curves are the band
  m <- matrix(sp$estimate__, 6, 5, byrow = TRUE)
  expect_equal(ce$x$lower__, apply(m, 2, stats::quantile, 0.025),
               ignore_attr = TRUE)
  ce2 <- conditional_effects(fit, "x:z", spaghetti = TRUE, band = "boot",
                             boot = 6, seed = 1, resolution = 5)
  expect_true(all(grepl("^[0-9]+_", levels(attr(ce2[[1]],
                                                "spaghetti")$sample__))))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(ce, ask = FALSE))
})

test_that("a nonlinear predictor gets a delta-method Wald band", {
  d <- ceo_data()
  fit <- frm(bf(y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1, nl = TRUE),
             family = gaussian(), data = d)
  ce <- conditional_effects(fit, "z", resolution = 20)
  df <- ce$z
  nd <- df
  expect_equal(df$estimate__,
               as.vector(frm_linpred(fit, newdata = nd, type = "response")))
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
  se <- sqrt(rowSums((lb$A %*% lb$V) * lb$A))
  # an absolute gap, bounded by a fraction of the standard error's scale
  expect_lt(max(abs(df$se__ - se)), 1e-8 * max(se))
  expect_true(all(df$lower__ < df$estimate__ & df$estimate__ < df$upper__))
  # the band that runs through one predictor is the only Wald band there
  expect_error(conditional_effects(fit, "z", band = "profile"),
               "nonlinear")
})

test_that("the nonlinear Wald band is refused where one Jacobian is short", {
  set.seed(4)
  d <- data.frame(z = stats::runif(200))
  d$y <- stats::rpois(200, 2 * exp(0.5 * d$z)) * stats::rbinom(200, 1, 0.8)
  fit <- frm(bf(y ~ a * exp(b * z), a ~ 1, b ~ 1, nl = TRUE),
             family = zero_inflated_poisson(), data = d)
  expect_error(conditional_effects(fit, "z"),
               "no Wald band for this display of a nonlinear predictor")
  # the predictor itself is one predictor
  ce <- conditional_effects(fit, "z", dpar = "mu", resolution = 10)
  expect_true(all(is.finite(ce$z$se__)))
})

test_that("a boot band at an observed group belongs to that group", {
  set.seed(21)
  d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:8, 10)))
  d$y <- stats::rnorm(80, 1 + 0.5 * d$x + stats::rnorm(8, 0, 2)[d$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  rg <- ranef(fit)$g
  re <- if (length(dim(rg)) == 3L) rg[, 1, 1] else rg[, 1]
  lv <- which.max(abs(re))
  at <- function(l, ...) {
    conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                        conditions = data.frame(g = factor(l, levels = 1:8)),
                        ...)
  }
  w <- at(lv)$x
  bo <- at(lv, band = "boot", boot = 40, seed = 5)
  b <- bo$x
  # the refits' curves center on the level's estimate, not on the
  # population curve a mode away; the base build was 0.86 modes off
  m <- colMeans(attr(bo, "boot")$t)
  expect_lt(max(abs(m - w$estimate__)), 0.3 * abs(re[lv]))
  expect_true(all(b$lower__ <= w$estimate__ & w$estimate__ <= b$upper__))
  # and two levels get two bands (the base build gave both level one's)
  other <- at(setdiff(1:8, lv)[1L], band = "boot", boot = 40, seed = 5)$x
  expect_gt(max(abs(other$lower__ - b$lower__)),
            0.1 * max(abs(other$estimate__ - w$estimate__)))
})
