# conditional_smooths() and posterior_average() on draws. The values
# are checked draw by draw against the fit machinery at that draw, and
# against brms at the same parameter vectors in
# dev/postfit2-brms-compare.R of the frmtmb repository.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

pf_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(11)
      dd <- data.frame(x = stats::runif(80), z = stats::rnorm(80))
      dd$y <- sin(2 * pi * dd$x) + 0.4 * dd$z + stats::rnorm(80, 0, 0.3)
      f1 <- frm(bf(y ~ z + s(x)), family = gaussian(), data = dd)
      f2 <- frm(bf(y ~ s(x)), family = gaussian(), data = dd)
      d1 <- suppressWarnings(suppressMessages(
        frm_sample(f1, chains = 1, iter = 300, refresh = 0, seed = 3)))
      d2 <- suppressWarnings(suppressMessages(
        frm_sample(f2, chains = 1, iter = 300, refresh = 0, seed = 4)))
      cache <<- list(dd = dd, f1 = f1, d1 = d1, d2 = d2)
    }
    cache
  }
})

test_that("conditional_smooths() on draws summarizes the per-draw terms", {
  cs <- pf_case()
  out <- conditional_smooths(cs$d1, resolution = 15, ndraws = 12,
                             spaghetti = TRUE)
  expect_identical(names(out), "mu: s(x)")
  fr <- out[[1]]
  sp <- attr(fr, "spaghetti")
  expect_equal(nrow(sp), 12 * 15)
  expect_identical(levels(sp$sample__), as.character(1:12))
  m <- matrix(sp$estimate__, 12, 15, byrow = TRUE)
  # brms's robust summary: median, MAD and quantiles of the curves
  expect_equal(fr$estimate__, apply(m, 2, stats::median))
  expect_equal(fr$se__, apply(m, 2, stats::mad))
  expect_equal(fr$upper__, apply(m, 2, stats::quantile, 0.975),
               ignore_attr = TRUE)
  # a curve is the fit's term at that draw: a contrast of its linear
  # predictor over a grid that varies x alone
  rows <- round(seq(1, ndraws(cs$d1), length.out = 12))
  idx <- frmtmb.sample:::draws_par_index(cs$d1$fit)
  for (k in c(1, 7)) {
    fk <- frmtmb.sample:::draws_fit_at(cs$d1, rows[k], idx)
    eta <- frm_linpred(fk, newdata = data.frame(x = fr$x, z = 0),
                       type = "link")
    # an absolute gap, bounded by a fraction of the curve's own scale
    expect_lt(max(abs((m[k, ] - m[k, 1]) - (eta - eta[1]))),
              1e-8 * max(abs(eta)))
  }
  # the grid is the fit method's
  fit_fr <- conditional_smooths(cs$f1, resolution = 15)[[1]]
  band <- c("estimate__", "se__", "lower__", "upper__")
  expect_identical(fr[setdiff(names(fr), band)],
                   fit_fr[setdiff(names(fit_fr), band)])
  expect_error(conditional_smooths(cs$d1, smooths = "s3"),
               "No valid smooth terms found in the model")
})

test_that("conditional_effects() on draws takes spaghetti and surface", {
  cs <- pf_case()
  ce <- conditional_effects(cs$d1, "x", resolution = 6, ndraws = 4,
                            spaghetti = TRUE)
  expect_equal(nrow(attr(ce$x, "spaghetti")), 4 * 6)
  expect_error(conditional_effects(cs$d1, "x", spaghetti = TRUE,
                                   surface = TRUE),
               "Cannot use 'spaghetti' and 'surface' at the same time")
  sf <- conditional_effects(cs$d1, "x:z", surface = TRUE, resolution = 5,
                            ndraws = 4)
  expect_true(attr(sf[[1]], "surface"))
  expect_equal(nrow(sf[[1]]), 25)
  far <- mgcv::exclude.too.far(sf[[1]]$x, sf[[1]]$z, cs$dd$x, cs$dd$z,
                               dist = 0.2)
  tf <- conditional_effects(cs$d1, "x:z", surface = TRUE, resolution = 5,
                            ndraws = 4, too_far = 0.2)
  expect_equal(nrow(tf[[1]]), sum(!far))
  sel <- conditional_effects(cs$d1, "x", ndraws = 4, resolution = 4,
                             select_points = 0.1)
  expect_lt(nrow(attr(sel$x, "points")), nrow(cs$dd))
})

test_that("posterior_average() takes draws in proportion to the weights", {
  cs <- pf_case()
  vars <- c("b_Intercept", "sigma")
  pa <- posterior_average(cs$d1, cs$d2, variable = vars,
                          weights = c(0.3, 0.7), seed = 5)
  n <- ndraws(cs$d1)
  expect_equal(dim(pa), c(n, 2))
  expect_identical(names(pa), vars)
  nd <- attr(pa, "ndraws")
  expect_equal(unname(nd), c(round(0.3 * n), n - round(0.3 * n)))
  expect_equal(unname(attr(pa, "weights")), c(0.3, 0.7))
  expect_identical(names(nd), c("cs$d1", "cs$d2"))
  # the rows are the chosen draws of each model, in model order
  set.seed(5)
  r1 <- sort(sample(seq_len(n), nd[1]))
  r2 <- sort(sample(seq_len(ndraws(cs$d2)), nd[2]))
  ref <- rbind(as.data.frame(cs$d1, variable = vars, draw = r1),
               as.data.frame(cs$d2, variable = vars, draw = r2))
  expect_equal(as.matrix(pa), as.matrix(ref), ignore_attr = TRUE)
  # a weight of zero takes nothing from that model
  pz <- posterior_average(cs$d1, cs$d2, variable = vars,
                          weights = c(1, 0), ndraws = 10, seed = 1)
  expect_equal(unname(attr(pz, "ndraws")), c(10, 0))
  # a variable only one model has needs `missing`
  expect_error(posterior_average(cs$d1, cs$d2, variable = "b_z",
                                 weights = c(1, 1)),
               "cannot be found in all of the models")
  pm <- posterior_average(cs$d1, cs$d2, variable = "b_z",
                          weights = c(1, 1), missing = 0, seed = 2)
  expect_equal(sum(pm$b_z == 0), unname(attr(pm, "ndraws")[2]))
  expect_error(posterior_average(cs$d1, cs$d2, weights = "kfold"),
               "weights = \"kfold\"")
  expect_error(posterior_average(cs$d1, cs$f1, weights = c(1, 1)),
               "averages draws objects")
})

test_that("posterior_average()'s stacking weights are loo's", {
  skip_if_not_installed("loo")
  cs <- pf_case()
  pa <- suppressWarnings(posterior_average(cs$d1, cs$d2,
                                           variable = "sigma", seed = 1))
  ref <- suppressWarnings(loo::loo_model_weights(
    list(loo(cs$d1), loo(cs$d2)), method = "stacking"))
  expect_equal(unname(attr(pa, "weights")),
               as.numeric(ref) / sum(as.numeric(ref)))
})

test_that("a curve at an observed group reads that group's draws", {
  # hand-built draws around the ML estimate: no sampler needed, and the
  # draws of one level can be moved by a known amount
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(1)
  M <- matrix(rep(est, each = 5) + stats::rnorm(5 * length(est), 0, 0.05),
              5, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  ds <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
  at <- function(o, l, ...) {
    conditional_effects(o, "x", resolution = 3, re_formula = NULL,
                        conditions = data.frame(g = factor(l, levels = 1:6)),
                        ...)$x
  }
  for (l in 1:2) {
    moved <- ds
    col <- sprintf("r_g[%d,Intercept]", l)
    moved$draws[, col] <- moved$draws[, col] + 10
    # the base build left the curve where it was, drawing a new group
    expect_equal(at(moved, l)$estimate__ - at(ds, l)$estimate__,
                 rep(10, 3))
    sp0 <- attr(at(ds, l, spaghetti = TRUE), "spaghetti")
    sp1 <- attr(at(moved, l, spaghetti = TRUE), "spaghetti")
    expect_equal(sp1$estimate__ - sp0$estimate__, rep(10, nrow(sp0)))
  }
  # levels one and two are different curves
  expect_false(isTRUE(all.equal(at(ds, 1)$estimate__, at(ds, 2)$estimate__)))
})

# hand-built draws around the ML estimate, as above: a known parameter
# vector per draw, so a level's effect can be moved by a known amount
hand_draws <- function(fit, n = 5, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}

test_that("on draws, an unseen or unset row is a new group, per row", {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- hand_draws(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd),
                   n = 40)
  at <- function(cond) {
    conditional_effects(ds, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, seed = 1)$x
  }
  new_g <- at(list())
  pop <- conditional_effects(ds, "x", resolution = 3)$x
  # the same draws of one new level, whatever the level is called; the
  # round-0 build gave "99" the population band, or refused it
  unseen <- at(list(g = "99"))
  expect_identical(as.character(unseen$g), rep("99", 3))
  expect_equal(unseen$lower__, new_g$lower__)
  expect_false(isTRUE(all.equal(unseen$lower__, pop$lower__)))
  cond <- data.frame(g = factor(c("2", NA), levels = 1:6))
  rownames(cond) <- c("lev2", "unset")
  mixed <- at(cond)
  expect_equal(mixed$lower__[mixed$cond__ == "unset"], new_g$lower__)
  expect_equal(mixed$lower__[mixed$cond__ == "lev2"],
               at(list(g = factor("2", levels = 1:6)))$lower__)
})

test_that("on draws, a nested term set at its parent reads the parent", {
  set.seed(25)
  d <- expand.grid(g = factor(1:8), h = factor(1:4), r = 1:5)
  d$x <- stats::rnorm(nrow(d))
  gh <- as.integer(interaction(d$g, d$h))
  d$y <- stats::rnorm(nrow(d), 1 + 0.5 * d$x +
                        stats::rnorm(8, 0, 2)[d$g] +
                        stats::rnorm(32, 0, 0.7)[gh])
  ds <- hand_draws(frm(bf(y ~ x + (1 | g / h)), family = gaussian(),
                       data = d), n = 40)
  at <- function(o, cond) {
    conditional_effects(o, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, seed = 1)$x
  }
  moved <- ds
  moved$draws[, "r_g[5,Intercept]"] <- moved$draws[, "r_g[5,Intercept]"] + 10
  cond <- list(g = factor("5", levels = 1:8))
  expect_equal(at(moved, cond)$estimate__ - at(ds, cond)$estimate__,
               rep(10, 3))
  # and the new g:h level's spread is in the band: wider than at an
  # observed g:h (the round-0 build gave the two the same width)
  both <- list(g = factor("5", levels = 1:8), h = factor("1", levels = 1:4))
  w_new <- at(ds, cond)
  w_obs <- at(ds, both)
  expect_gt(min((w_new$upper__ - w_new$lower__) /
                  (w_obs$upper__ - w_obs$lower__)), 1.5)
})

test_that("on draws, a grouping variable varied as an effect is observed", {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- hand_draws(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
  xg <- conditional_effects(ds, "x:g", resolution = 3,
                            re_formula = NULL)[["x:g"]]
  one <- conditional_effects(ds, "x", resolution = 3, re_formula = NULL,
                             conditions = list(g = factor("4", levels = 1:6)))
  rows <- as.character(xg$g) == "4"
  expect_equal(xg$estimate__[rows], one$x$estimate__)
  expect_equal(xg$upper__[rows], one$x$upper__)
})

test_that("on draws, a gr(g, by = f) row reads its own f level's term", {
  set.seed(44)
  d <- data.frame(x = stats::rnorm(240), g = factor(rep(1:12, 20)))
  d$f <- factor(ifelse(as.integer(d$g) <= 6, "a", "b"))
  sdg <- ifelse(1:12 <= 6, 0.3, 1.5)
  d$y <- stats::rnorm(240, 1 + 0.5 * d$x + stats::rnorm(12, 0, sdg)[d$g],
                      0.5)
  ds <- hand_draws(frm(bf(y ~ x + f + (1 | gr(g, by = f))),
                       family = gaussian(), data = d), n = 40)
  at <- function(o, cond) {
    conditional_effects(o, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, seed = 1)$x
  }
  # the round-1 build refused both calls, and a draws display has no
  # other band
  moved <- ds
  moved$draws[, "r_g[3,Intercept]"] <- moved$draws[, "r_g[3,Intercept]"] + 10
  obs <- list(f = "a", g = "3")
  expect_equal(at(moved, obs)$estimate__ - at(ds, obs)$estimate__,
               rep(10, 3))
  # a new g within f = a is drawn with f = a's small sd
  new_a <- at(ds, list(f = "a"))
  new_b <- at(ds, list(f = "b"))
  expect_lt(max((new_a$upper__ - new_a$lower__) /
                  (new_b$upper__ - new_b$lower__)), 0.5)
})

test_that("on draws, an mm() term reads a level per member", {
  set.seed(43)
  d <- data.frame(x = stats::rnorm(200),
                  g1 = factor(sample(1:10, 200, TRUE)),
                  g2 = factor(sample(1:10, 200, TRUE)))
  u <- stats::rnorm(10, 0, 1)
  d$y <- stats::rnorm(200, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]), 0.5)
  ds <- hand_draws(frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(),
                       data = d), n = 40)
  at <- function(o, cond) {
    conditional_effects(o, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, seed = 1)$x
  }
  moved <- ds
  col <- "r_mmg1g2[2,Intercept]"
  moved$draws[, col] <- moved$draws[, col] + 10
  # both members at level 2: weight 1/2 each, so the full shift
  both <- list(g1 = "2", g2 = "2")
  expect_equal(at(moved, both)$estimate__ - at(ds, both)$estimate__,
               rep(10, 3))
  # one member at level 2 and one unset: half the shift, and the unset
  # member's new level is in the band
  one <- list(g1 = "2")
  expect_equal(at(moved, one)$estimate__ - at(ds, one)$estimate__,
               rep(5, 3))
  # nothing set: a new group, not the first observed one (round 1)
  new <- at(ds, list())
  expect_true(all(is.na(new$g1) & is.na(new$g2)))
  expect_gt(min((new$upper__ - new$lower__) /
                  (at(ds, both)$upper__ - at(ds, both)$lower__)), 1.5)
  # the two unset members are ONE new level with their weights added,
  # the Wald band's reading and brms's: over 200 draws the band is
  # 0.94 to 0.97 of the Wald band, and 0.60 to 0.66 with the members
  # drawn as two levels (dev/postfit2-rev-p2-mmsplit.R)
  d200 <- at(hand_draws(ds$fit, n = 200), list())
  wald <- conditional_effects(ds$fit, "x", resolution = 3,
                              re_formula = NULL)$x
  r <- (d200$upper__ - d200$lower__) / (wald$upper__ - wald$lower__)
  expect_gt(min(r), 0.8)
  expect_lt(max(r), 1.25)
})

test_that("on draws, mm() beside (1 | g1) draws both new levels", {
  # refused before (review P2-M1): the plain term placed g1 first
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
  ds <- hand_draws(fit, n = 200)
  for (cond in list(list(), list(g1 = "99"))) {
    a <- conditional_effects(ds, "x", resolution = 3, re_formula = NULL,
                             conditions = cond, seed = 1)$x
    w <- conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                             conditions = cond)$x
    r <- (a$upper__ - a$lower__) / (w$upper__ - w$lower__)
    expect_gt(min(r), 0.8)
    expect_lt(max(r), 1.25)
  }
})

test_that("on draws, a new level is drawn once and shared by every panel", {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60),
                   g = factor(rep(1:6, 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x + 0.3 * dd$z +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- hand_draws(frm(bf(y ~ x + z + (1 | g)), family = gaussian(),
                       data = dd), n = 40)
  # both panels hold one row at the same point and read the same new g
  ce <- conditional_effects(ds, c("x", "z"), re_formula = NULL,
                            int_conditions = list(x = mean(dd$x),
                                                  z = mean(dd$z)),
                            seed = 1)
  expect_identical(ce$x$estimate__, ce$z$estimate__)
  expect_identical(ce$x$lower__, ce$z$lower__)
  expect_identical(ce$x$upper__, ce$z$upper__)
})

test_that("on draws, crossed terms at an unseen combination answer", {
  # (1 | g) + (1 | h) + (1 | g:h) with g and h at observed levels the
  # data never has together. frmtmb.sample on 0.66.0 refused it ("cannot
  # draw a new level of the group-level term (1 | g:h)"); against brms
  # 2.23.0 at the same draws and seed the band agrees to 4.4e-16
  # (dev/ceplot-log/crossed-brms.txt)
  set.seed(49)
  dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
  dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
  dc$x <- stats::rnorm(nrow(dc))
  dc$y <- stats::rnorm(nrow(dc), 1 + 0.5 * dc$x + stats::rnorm(6)[dc$g] +
                         stats::rnorm(5)[dc$h] +
                         stats::rnorm(30, 0, 0.7)[
                           as.integer(interaction(dc$g, dc$h))], 0.5)
  fc <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
            data = dc)
  ds <- hand_draws(fc, n = 200)
  at <- function(o, cond, ...) {
    conditional_effects(o, "x", resolution = 3, re_formula = NULL,
                        conditions = cond, ...)$x
  }
  new <- list(g = "1", h = "1")
  a <- at(ds, new, seed = 1)
  w <- at(fc, new)
  expect_identical(as.character(a$g), rep("1", 3))
  r <- (a$upper__ - a$lower__) / (w$upper__ - w$lower__)
  expect_gt(min(r), 0.7)
  expect_lt(max(r), 1.4)
  # the g:h level is drawn: an observed combination is much narrower
  o <- at(ds, list(g = "1", h = "2"), seed = 1)
  expect_lt(max((o$upper__ - o$lower__) / (a$upper__ - a$lower__)), 0.5)
  # moving the observed g's effect moves the curve by that much: g is
  # read at its own level, not moved to a placeholder
  moved <- ds
  moved$draws[, "r_g[1,Intercept]"] <- moved$draws[, "r_g[1,Intercept]"] + 10
  expect_equal(at(moved, new, seed = 1)$estimate__ - a$estimate__,
               rep(10, 3))
})

test_that("on draws, an mm() term with a by variable draws new members", {
  # frmtmb.sample on 0.66.0 refused every new member of a by-split
  # mm() term; each is now drawn in its own by-level's block
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
  ds <- hand_draws(fit, n = 200)
  for (cond in list(list(f1 = "a", f2 = "b"), list(f1 = "a", f2 = "a"),
                    list(g1 = "2", f1 = "a", f2 = "b"))) {
    a <- conditional_effects(ds, "x", resolution = 3, re_formula = NULL,
                             conditions = cond, seed = 1)$x
    w <- conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                             conditions = cond)$x
    r <- (a$upper__ - a$lower__) / (w$upper__ - w$lower__)
    expect_gt(min(r), 0.7)
    expect_lt(max(r), 1.4)
  }
})
