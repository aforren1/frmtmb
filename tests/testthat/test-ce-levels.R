# Which group a conditional_effects() grid row reads on a FIT: brms's
# rule, per row and per group-level term (allow_new_levels = TRUE in
# brms's display), under the Wald band and the bootstrap band. The draws
# method's side is frmtmb.sample's test-postfit-draws.R.

lvl_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(21)
      d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:8, 10)))
      d$y <- stats::rnorm(80, 1 + 0.5 * d$x + stats::rnorm(8, 0, 2)[d$g], 1)
      cache <<- list(d = d, fit = frm(bf(y ~ x + (1 | g)),
                                      family = gaussian(), data = d))
    }
    cache
  }
})
width <- function(df) df$upper__ - df$lower__
at_g <- function(fit, cond, ...) {
  conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                      conditions = cond, ...)$x
}

test_that("a level the fit never saw is a new group, and keeps its name", {
  s <- lvl_fit()
  new_g <- at_g(s$fit, list())
  unseen <- at_g(s$fit, list(g = "99"))
  expect_identical(as.character(unseen$g), rep("99", 3))
  expect_equal(unseen$lower__, new_g$lower__)
  expect_equal(unseen$upper__, new_g$upper__)
  # the bootstrap draws that new level too: the round-0 build gave the
  # POPULATION band here, 0.17 of the Wald width
  b <- at_g(s$fit, list(g = "99"), band = "boot", boot = 30, seed = 5)
  expect_gt(min(width(b) / width(unseen)), 0.5)
})

test_that("rows mixing a level and NA read one each, or boot says why not", {
  s <- lvl_fit()
  cond <- data.frame(g = factor(c("2", NA), levels = levels(s$d$g)))
  rownames(cond) <- c("lev2", "unset")
  mixed <- at_g(s$fit, cond)
  expect_equal(mixed$upper__[mixed$cond__ == "unset"],
               at_g(s$fit, list())$upper__)
  expect_equal(mixed$upper__[mixed$cond__ == "lev2"],
               at_g(s$fit, list(g = factor("2", levels = 1:8)))$upper__)
  expect_error(at_g(s$fit, cond, band = "boot", boot = 5, seed = 1),
               "observed level in some grid rows and at a new level")
})

test_that("a nested term set at its parent draws a new level within it", {
  set.seed(25)
  d <- expand.grid(g = factor(1:8), h = factor(1:4), r = 1:5)
  d$x <- stats::rnorm(nrow(d))
  gh <- as.integer(interaction(d$g, d$h))
  d$y <- stats::rnorm(nrow(d), 1 + 0.5 * d$x +
                        stats::rnorm(8, 0, 2)[d$g] +
                        stats::rnorm(32, 0, 0.7)[gh])
  fit <- frm(bf(y ~ x + (1 | g / h)), family = gaussian(), data = d)
  rg <- ranef(fit)$g
  re <- if (length(dim(rg)) == 3L) rg[, 1, 1] else rg[, 1]
  cond <- list(g = factor(which.max(abs(re)), levels = 1:8))
  w <- at_g(fit, cond)
  b <- at_g(fit, cond, band = "boot", boot = 40, seed = 5)
  # the round-0 build lost the new g:h variance: 0.25 of the Wald width
  r <- width(b) / width(w)
  expect_gt(min(r), 0.6)
  expect_lt(max(r), 1.6)
  # and the estimate is g = 2's: its own mode is read
  pop <- conditional_effects(fit, "x", resolution = 3)$x
  expect_gt(max(abs(w$estimate__ - pop$estimate__)),
            0.1 * max(abs(width(w))))
})

test_that("crossed terms, one set: the boot band sits on its estimate", {
  set.seed(22)
  d <- data.frame(x = stats::rnorm(160), g = factor(rep(1:8, 20)),
                  h = factor(rep(1:10, each = 16)))
  d$y <- stats::rnorm(160, 1 + 0.5 * d$x + stats::rnorm(8, 0, 2)[d$g] +
                        stats::rnorm(10, 0, 0.5)[d$h], 1)
  fit <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = d)
  rg <- ranef(fit)$g
  re <- if (length(dim(rg)) == 3L) rg[, 1, 1] else rg[, 1]
  lv <- which.max(abs(re))
  cond <- list(g = factor(lv, levels = 1:8))
  w <- at_g(fit, cond)
  bo <- conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                            conditions = cond, band = "boot", boot = 40,
                            seed = 5)
  # the round-0 build: 1.035 modes off, and the band excluded its own
  # estimate at 2 of 3 points
  m <- colMeans(attr(bo, "boot")$t)
  expect_lt(max(abs(m - w$estimate__)), 0.4 * abs(re[lv]))
  expect_true(all(bo$x$lower__ <= w$estimate__ &
                    w$estimate__ <= bo$x$upper__))
})

test_that("re_formula changes nothing on a model with no group term", {
  set.seed(24)
  d <- data.frame(x = stats::runif(100))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(100, 0, 0.3)
  fit <- frm(bf(y ~ s(x)), family = gaussian(), data = d)
  a <- conditional_effects(fit, "x", resolution = 5, band = "boot",
                           boot = 20, seed = 5)$x
  b <- conditional_effects(fit, "x", resolution = 5, re_formula = NULL,
                           band = "boot", boot = 20, seed = 5)$x
  # the round-0 build shrank the NULL band to 0.06 to 0.09 of this one
  expect_identical(b$lower__, a$lower__)
  expect_identical(b$upper__, a$upper__)
})

test_that("a group only another parameter reads does not refuse the band", {
  set.seed(7)
  d <- data.frame(z = stats::runif(120), g = factor(rep(1:6, 20)))
  d$y <- 2 * exp(0.8 * d$z) +
    stats::rnorm(120, 0, exp(-1 + stats::rnorm(6, 0, 0.3)[d$g]))
  fit <- frm(bf(y ~ a * exp(b * z), a ~ 1, b ~ 1, sigma ~ (1 | g),
                nl = TRUE), family = gaussian(), data = d)
  ce <- conditional_effects(fit, "z", resolution = 4, re_formula = NULL)
  expect_true(all(is.finite(ce$z$se__)))
  # the guard still fires where the displayed body reads the new group
  fg <- frm(bf(y ~ a * exp(b * z), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
            family = gaussian(), data = d)
  expect_error(conditional_effects(fg, "z", resolution = 4,
                                   re_formula = NULL),
               "no Wald band for a nonlinear predictor at a new group")
})

test_that("a population bootstrap is not reused for an observed group", {
  s <- lvl_fit()
  # the same grid and the same conditions, read at the population level
  # and at an observed group: only the simulation tells them apart, so
  # only the key's `kept` field can refuse the swap
  cond <- list(g = factor("1", levels = levels(s$d$g)))
  pop <- conditional_effects(s$fit, "x", resolution = 3, conditions = cond,
                             band = "boot", boot = 5, seed = 2)
  expect_identical(pop$x$g, at_g(s$fit, cond)$g)
  expect_error(at_g(s$fit, cond, band = "boot", boot = attr(pop, "boot")),
               "refits of a different simulation")
  # and the other way round
  obs <- conditional_effects(s$fit, "x", resolution = 3, re_formula = NULL,
                             conditions = cond, band = "boot", boot = 5,
                             seed = 2)
  expect_error(conditional_effects(s$fit, "x", resolution = 3,
                                   conditions = cond, band = "boot",
                                   boot = attr(obs, "boot")),
               "refits of a different simulation")
})

test_that("a grouping variable varied as an effect is an observed group", {
  s <- lvl_fit()
  xg <- conditional_effects(s$fit, "x:g", resolution = 3, re_formula = NULL,
                            band = "boot", boot = 20, seed = 3)[["x:g"]]
  one <- at_g(s$fit, list(g = factor("3", levels = 1:8)), band = "boot",
              boot = 20, seed = 3)
  # the same simulated data and refits, so the level-3 rows are the
  # level-3 call's band; the round-0 build drew a new level for them
  rows <- as.character(xg$g) == "3"
  expect_equal(xg$lower__[rows], one$lower__)
  expect_equal(xg$upper__[rows], one$upper__)
})

test_that("the bootstrap band holds the smooths at their fit", {
  # user decision, 2026-09-29: the display's bootstrap never redraws a
  # smooth from its prior; the base build's band on y ~ s(x) was 10 to
  # 18 times the Wald band, under both re_formula values
  in_band <- function(fit, ...) {
    w <- conditional_effects(fit, "x", resolution = 5, ...)$x
    b <- conditional_effects(fit, "x", resolution = 5, band = "boot",
                             boot = 40, seed = 5, ...)$x
    r <- (b$upper__ - b$lower__) / (w$upper__ - w$lower__)
    expect_gt(min(r), 0.5)
    expect_lt(max(r), 2)
    expect_true(all(b$lower__ <= w$estimate__ &
                      w$estimate__ <= b$upper__))
  }
  set.seed(24)
  d <- data.frame(x = stats::runif(100))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(100, 0, 0.3)
  f1 <- frm(bf(y ~ s(x)), family = gaussian(), data = d)
  in_band(f1)
  in_band(f1, re_formula = NULL)
  # beside a group-level term: at the population, at a new group and
  # at an observed one
  set.seed(26)
  d2 <- data.frame(x = stats::runif(160), g = factor(rep(1:8, 20)))
  d2$y <- sin(2 * pi * d2$x) + stats::rnorm(8)[d2$g] +
    stats::rnorm(160, 0, 0.3)
  f2 <- frm(bf(y ~ s(x) + (1 | g)), family = gaussian(), data = d2)
  in_band(f2)
  in_band(f2, re_formula = NULL)
  in_band(f2, re_formula = NULL,
          conditions = list(g = factor("3", levels = 1:8)))
})

test_that("a gr(g, by = f) row reads only the term of its own f level", {
  set.seed(44)
  d <- data.frame(x = stats::rnorm(240), g = factor(rep(1:12, 20)))
  d$f <- factor(ifelse(as.integer(d$g) <= 6, "a", "b"))
  sdg <- ifelse(1:12 <= 6, 0.3, 1.5)
  d$y <- stats::rnorm(240, 1 + 0.5 * d$x + stats::rnorm(12, 0, sdg)[d$g],
                      0.5)
  fit <- frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(),
             data = d)
  ratio <- function(cond) {
    w <- at_g(fit, cond)
    b <- at_g(fit, cond, band = "boot", boot = 40, seed = 5)
    expect_true(all(b$lower__ <= w$estimate__ & w$estimate__ <= b$upper__))
    width(b) / width(w)
  }
  # the round-1 build refused both, reading the f = b term in every row
  r_obs <- ratio(list(f = "a", g = "3"))
  expect_gt(min(r_obs), 0.5)
  expect_lt(max(r_obs), 2)
  r_new <- ratio(list(f = "a"))
  expect_gt(min(r_new), 0.5)
  expect_lt(max(r_new), 2)
  # a new g within f = a is drawn with f = a's small sd, not f = b's
  new_a <- at_g(fit, list(f = "a"), band = "boot", boot = 40, seed = 5)
  new_b <- at_g(fit, list(f = "b"), band = "boot", boot = 40, seed = 5)
  expect_lt(max(width(new_a) / width(new_b)), 0.5)
})

test_that("a new level's placeholder never moves a fixed effect's column", {
  # y ~ x + trt + (1 | trt:subj) with trt set and subj unset: the new
  # trt:subj level has to be placed at trt = "b". With nothing locked,
  # it moved trt to "a" and the curve with it (boot mean 0.288, not
  # 1.092); dev/postfit2-rev-p1-nolock.R
  set.seed(41)
  d <- expand.grid(rep = 1:4, trt = factor(c("a", "b")), subj = factor(1:12))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d), 1 + 0.5 * d$x + (d$trt == "b") +
                        stats::rnorm(24, 0, 1)[interaction(d$trt, d$subj)],
                      0.5)
  fit <- frm(bf(y ~ x + trt + (1 | trt:subj)), family = gaussian(), data = d)
  at_trt <- function(l, ...) {
    conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                        conditions = list(trt = l), ...)
  }
  b <- at_trt("b", band = "boot", boot = 40, seed = 5)
  m <- colMeans(attr(b, "boot")$t, na.rm = TRUE)
  est_b <- at_trt("b")$x$estimate__
  est_a <- at_trt("a")$x$estimate__
  expect_true(all(abs(m - est_b) < abs(m - est_a)))
  expect_identical(as.character(b$x$trt), rep("b", 3))
})

test_that("an mm() term's members are grouping variables", {
  set.seed(43)
  d <- data.frame(x = stats::rnorm(200),
                  g1 = factor(sample(1:10, 200, TRUE)),
                  g2 = factor(sample(1:10, 200, TRUE)))
  u <- stats::rnorm(10, 0, 1)
  d$y <- stats::rnorm(200, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]), 0.5)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(), data = d)
  mm_at <- function(cond, ...) {
    conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, ...)$x
  }
  pop <- conditional_effects(fit, "x", resolution = 3)$x
  # nothing set: a new group, not the first observed one silently
  new <- mm_at(list())
  expect_true(all(is.na(new$g1) & is.na(new$g2)))
  expect_equal(new$estimate__, pop$estimate__)
  expect_true(all(width(new) > width(pop)))
  # an observed group: the boot band sits on its estimate and has the
  # Wald band's width (the round-1 build redrew the term: 5 to 6 times
  # as wide, near the population curve)
  both2 <- list(g1 = "2", g2 = "2")
  w <- mm_at(both2)
  b <- mm_at(both2, band = "boot", boot = 60, seed = 5)
  expect_true(all(b$lower__ <= w$estimate__ & w$estimate__ <= b$upper__))
  expect_lt(max(width(b) / width(w)), 2)
  # a new group on boot: two unset members are one new level with their
  # weights added, as the Wald band reads them. Drawn as two levels of
  # weight 1/2 each, the band was 25 to 38% narrower
  # (dev/postfit2-rev-p2-mmsplit.R): boot over Wald 0.67 to 0.69
  bn <- mm_at(list(), band = "boot", boot = 60, seed = 5)
  expect_gt(min(width(bn) / width(new)), 0.8)
  expect_lt(max(width(bn) / width(new)), 1.25)
  # one member observed and one new in the same row
  expect_error(mm_at(list(g1 = "2"), band = "boot", boot = 5, seed = 1),
               "different members of one row")
})

test_that("mm() beside (1 | g1) shares g1's placeholder when both are new", {
  # the plain term placed g1 first, and the mm() term was then refused
  # (review P2-M1); the two terms' coefficients are separate slots
  set.seed(61)
  d <- data.frame(x = stats::rnorm(240),
                  g1 = factor(sample(1:10, 240, TRUE)),
                  g2 = factor(sample(1:10, 240, TRUE)))
  u <- stats::rnorm(10, 0, 1)
  v <- stats::rnorm(10, 0, 0.7)
  d$y <- stats::rnorm(240, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]) +
                        v[d$g1], 0.5)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2)) + (1 | g1)), family = gaussian(),
             data = d)
  for (cond in list(list(), list(g1 = "99"))) {
    w <- at_g(fit, cond)
    b <- at_g(fit, cond, band = "boot", boot = 40, seed = 5)
    expect_gt(min(width(b) / width(w)), 0.5)
    expect_lt(max(width(b) / width(w)), 2)
  }
})

test_that("an mm() term with a by variable holds its observed members", {
  set.seed(45)
  d <- data.frame(x = stats::rnorm(300),
                  g1 = factor(sample(1:10, 300, TRUE)),
                  g2 = factor(sample(1:10, 300, TRUE)))
  fl <- rep(c("a", "b"), each = 5)
  d$f1 <- factor(fl[d$g1])
  d$f2 <- factor(fl[d$g2])
  u <- stats::rnorm(10, 0, 1)
  d$y <- stats::rnorm(300, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]), 0.5)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))),
             family = gaussian(), data = d)
  # members at levels of the two by-levels' blocks: both observed
  cond <- list(g1 = "2", g2 = "7", f1 = "a", f2 = "b")
  w <- at_g(fit, cond)
  b <- at_g(fit, cond, band = "boot", boot = 40, seed = 5)
  expect_true(all(b$lower__ <= w$estimate__ & w$estimate__ <= b$upper__))
  expect_lt(max(width(b) / width(w)), 2)
  expect_error(at_g(fit, list(f1 = "a", f2 = "b"), band = "boot", boot = 5,
                    seed = 1),
               "multi-membership term .* with a by variable")
})

test_that("a bootstrap is not reused under another re_formula", {
  set.seed(48)
  d <- data.frame(x = stats::rnorm(240), g = factor(rep(1:12, 20)))
  b0 <- stats::rnorm(12, 0, 1)
  b1 <- stats::rnorm(12, 0, 1)
  d$y <- stats::rnorm(240, 1 + b0[d$g] + (0.5 + b1[d$g]) * d$x, 0.5)
  fit <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)
  cond <- list(g = "3")
  a <- conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                           conditions = cond, band = "boot", boot = 5,
                           seed = 1)
  # the same grid and the same term held: only re_formula differs, and
  # ~ (1 | g) drops the slope the NULL draws carry
  expect_error(conditional_effects(fit, "x", resolution = 3,
                                   re_formula = ~ (1 | g),
                                   conditions = cond, band = "boot",
                                   boot = attr(a, "boot")),
               "different re_formula")
  again <- conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                               conditions = cond, band = "boot",
                               boot = attr(a, "boot"))
  expect_equal(again$x$upper__, a$x$upper__)
})

test_that("a smooth's factor at a level it has no curve for is refused", {
  set.seed(47)
  d <- data.frame(x = stats::runif(200), g = factor(rep(1:5, 40)))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(5, 0, 0.5)[d$g] +
    stats::rnorm(200, 0, 0.3)
  fs <- frm(bf(y ~ s(x, g, bs = "fs", k = 5)), family = gaussian(),
            data = d)
  # the round-1 build returned an NA estimate here without a word
  for (rf in list(NULL, NA)) {
    expect_error(conditional_effects(fs, "x", resolution = 3,
                                     re_formula = rf,
                                     conditions = list(g = "99")),
                 "sets g to \"99\", a level the fit never saw")
  }
  expect_true(all(is.finite(at_g(fs, list(g = "2"))$estimate__)))
  # beside (1 | g), re_formula = NULL leaves g unset
  fg <- frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | g)),
            family = gaussian(), data = d)
  expect_error(at_g(fg, list()), "leaves g unset (NA)", fixed = TRUE)
})

test_that("a new level is drawn once and shared by every panel", {
  s <- lvl_fit()
  set.seed(3)
  d <- s$d
  d$z <- stats::rnorm(nrow(d))
  fit <- frm(bf(y ~ x + z + (1 | g)), family = gaussian(), data = d)
  # both panels hold one row at the same point, x = mean(x) and
  # z = mean(z), and both read the same new level of g
  ce <- conditional_effects(fit, c("x", "z"), re_formula = NULL,
                            int_conditions = list(x = mean(d$x),
                                                  z = mean(d$z)),
                            band = "boot", boot = 20, seed = 4)
  t <- attr(ce, "boot")$t
  expect_identical(ncol(t), 2L)
  expect_identical(t[, 1], t[, 2])
  expect_identical(ce$x$upper__, ce$z$upper__)
})
