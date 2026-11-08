# A me(gr = ) level the fit did not see, on draws, under allow_new_levels:
# each draw takes the level's latent value from normal(x, sdx), as brms does
# and as fitted() integrates. Through 0.26.0 it was refused ("cannot predict a
# grouping level the fit did not see") where fitted() and brms answered: 8.051
# and 7.979 at one row, seed 7511. A group-level effect's new level stays
# refused. The draws here all sit at the estimate, so their spread is the
# latent value's alone.

mg_case <- function() {
  set.seed(7511)
  nid <- 40
  per <- 5
  id <- rep(seq_len(nid), each = per)
  lat <- stats::rnorm(nid)
  d <- data.frame(id = factor(id), g = factor(rep(1:8, length.out = nid)[id]),
                  x = stats::rnorm(nid * per))
  d$zo <- (lat + stats::rnorm(nid, 0, 0.4))[id]
  d$sz <- 0.4
  d$y <- stats::rpois(nrow(d), exp(0.5 + 0.6 * lat[id] + 0.3 * d$x +
                                     stats::rnorm(8, 0, 0.3)[d$g]))
  fit <- frm(bf(y ~ x + me(zo, sz, gr = id) + (1 | g)), family = poisson(),
             data = d)
  idx <- draws_par_index(fit)
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  row <- numeric(length(lab))
  for (cp in names(idx)) row[idx[[cp]]] <- fit$estimates[[cp]]
  m <- matrix(rep(row, each = 4000), 4000, dimnames = list(NULL, lab))
  ds <- structure(list(stanfit = NULL, draws = draws_to_natural(m, fit),
                       fit = fit), class = "frmtmb_draws")
  list(fit = fit, ds = ds)
}

test_that("a new me(gr = ) level on draws reads its draw, as fitted() does", {
  mc <- mg_case()
  nd <- data.frame(id = "new", g = "1", x = 0, zo = 0.3, sz = 0.4)
  set.seed(1)
  e <- posterior_epred(mc$ds, newdata = nd, allow_new_levels = TRUE)
  f <- fitted(mc$fit, newdata = nd, allow_new_levels = TRUE)
  mm <- frmtmb:::fitted_me_moments(mc$fit, nd, NULL, NULL, NULL, TRUE)
  se <- stats::sd(e[, 1]) / sqrt(nrow(e))
  expect_lt(abs(mean(e[, 1]) - f[1, "Estimate"]) / se, 4)
  expect_lt(abs(stats::sd(e[, 1]) / sqrt(mm$var) - 1), 4 / sqrt(2 * nrow(e)))
  # one level read by two rows is one draw
  nd2 <- rbind(nd, transform(nd, x = 1))
  e2 <- posterior_epred(mc$ds, newdata = nd2, allow_new_levels = TRUE)
  expect_gt(stats::cor(log(e2[, 1]), log(e2[, 2])), 0.999)
  # without the flag it is refused, in brms's words
  expect_error(posterior_epred(mc$ds, newdata = nd), "allow_new_levels")
  # a new level of the group-level term stays refused under the flag
  expect_error(posterior_epred(mc$ds, newdata = transform(nd, g = "new"),
                               allow_new_levels = TRUE),
               "cannot predict a grouping level")
})
