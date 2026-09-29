# A refit inside the package carries the FITTED model's ordinal
# threshold count. Recounting it from a leave-one-out subset or from a
# simulated replicate gives a parameter vector of a different length and
# a different meaning, silently. dev/thresrefit-findings.md has the
# before-and-after numbers and the scripts.

# n = 50, four categories, EXACTLY ONE row in the top category, so
# deleting that row removes the category from the data.
refit_one_top <- function(seed = 501, n = 50) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x),
              stats::plogis(2.3 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  top <- which(y == 4L)
  if (!length(top)) stop("the seed reached no top category")
  y[top[-1L]] <- 3L
  data.frame(x = x, y = y)
}

# 96 rows over three levels of g, level "a" alone reaching category 4
# and there on exactly one row.
refit_grouped_top <- function(seed = 502, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  ia <- which(g == "a" & y == 4L)
  if (!length(ia)) stop("the seed reached no top category in level a")
  y[ia[-1L]] <- 3L
  data.frame(x = x, g = g, y = y)
}

test_that("a leave-one-out refit keeps the fitted threshold count", {
  skip_on_cran()
  dd <- refit_one_top()
  itop <- which(dd$y == 4L)
  expect_length(itop, 1L)

  for (fam in list(cumulative(), sratio())) {
    fit <- suppressWarnings(frm(bf(y ~ x), family = fam, data = dd))
    expect_length(fit$estimates[["tau_raw"]], 3L)
    inf <- suppressWarnings(influence(fit, force = TRUE))
    # the whole table, not only the deleted row: a shorter refit used to
    # leave an NA behind and shift the rest of the row
    expect_false(anyNA(inf$fixed))
    expect_false(anyNA(cooks.distance(inf)))
    expect_identical(colnames(inf$fixed),
                     c("x", "tau_raw_1", "tau_raw_2", "tau_raw_3"))

    # the reference: the same data with the count PINNED by hand, read
    # through the accessor the influence table is built from, so the
    # order and the names are the same on both sides
    ref <- suppressWarnings(
      frm(bf(y | thres(3) ~ x), family = fam,
          data = dd[-itop, , drop = FALSE]))
    want <- frmtmb:::get_coef.frmtmb_fit(ref)
    got <- inf$fixed[itop, ]
    expect_identical(names(want), colnames(inf$fixed))
    # The top threshold bounds a category nobody chose once the row is
    # gone, so it is not identified: the two optimizations stop at
    # different points on a flat ridge. Every identified coefficient is
    # compared as a ratio to the spread the deletions themselves produce
    # in that column, which the run measures.
    ident <- setdiff(names(want), "tau_raw_3")
    spread <- apply(inf$fixed, 2, stats::sd)
    expect_lt(max(abs(got[ident] - want[ident]) / spread[ident]), 1e-3)
    # and the unidentified one is large on both sides, in units of the
    # largest identified threshold
    scale_id <- max(abs(want[c("tau_raw_1", "tau_raw_2")]))
    expect_gt(abs(got[["tau_raw_3"]]) / scale_id, 3)
    expect_gt(abs(want[["tau_raw_3"]]) / scale_id, 3)
  }
})

test_that("a leave-one-out refit keeps the per-level threshold counts", {
  skip_on_cran()
  dd <- refit_grouped_top()
  ia <- which(dd$g == "a" & dd$y == 4L)
  expect_length(ia, 1L)
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                              family = cumulative(), data = dd))
  th <- fit$spec$responses$y$family[["thres"]]
  expect_identical(th[["nthres"]], c(3L, 2L, 2L))

  inf <- suppressWarnings(influence(fit, force = TRUE))
  expect_false(anyNA(inf$fixed))
  expect_false(anyNA(cooks.distance(inf)))
  expect_identical(ncol(inf$fixed), 8L)

  dsub <- dd[-ia, , drop = FALSE]
  dsub$k <- as.integer(th[["nthres"]][match(as.character(dsub$g),
                                            th[["levels"]])])
  ref <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                              family = cumulative(), data = dsub))
  want <- frmtmb:::get_coef.frmtmb_fit(ref)
  got <- inf$fixed[ia, ]
  expect_identical(names(want), colnames(inf$fixed))
  # level a's third threshold is position 3 of the merged vector, and
  # it is the unidentified one; every other coefficient has to line up
  # POSITION BY POSITION, which is what the recount used to break by
  # shifting levels b and c down one slot
  unident <- paste0("tau_raw_", th[["nthres"]][1L])
  ident <- setdiff(names(want), unident)
  spread <- apply(inf$fixed, 2, stats::sd)
  expect_lt(max(abs(got[ident] - want[ident]) / spread[ident]), 1e-3)
  scale_id <- max(abs(want[setdiff(ident, "x")]))
  expect_gt(abs(got[[unident]]) / scale_id, 3)
})

test_that("thres_pin_of_fit() finds nothing to pin off an ordinal model", {
  skip_on_cran()
  set.seed(505)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- stats::rnorm(40, dd$x)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_null(frmtmb:::thres_pin_of_fit(fit))
  # and a pin that names no response leaves the values alone
  av <- list(weights = rep(1, 3))
  expect_identical(
    frmtmb:::thres_pin_apply(list(other = list(grouped = FALSE, nthres = 3)),
                             list(resp_name = "y"), av),
    av)
  # a count the caller wrote wins over the pinned COUNT
  av2 <- list(thres = 5)
  expect_identical(
    frmtmb:::thres_pin_apply(list(y = list(grouped = FALSE, nthres = 3)),
                             list(resp_name = "y"), av2)[["thres"]], 5)
  # but not over the RECODING, which is about what the codes mean. A
  # subset missing the interior category "b" codes its rows 1,2,3 for
  # a,c,d, and on the fitted scale those are 1,3,4.
  p <- list(y = list(grouped = FALSE, nthres = 3,
                     levels_y = c("a", "b", "c", "d")))
  rc <- frmtmb:::thres_pin_recode(p, list(resp_name = "y"), c(1, 2, 3),
                                  c("a", "c", "d"))
  expect_identical(rc$y, c(1L, 3L, 4L))
  expect_identical(rc$levels, c("a", "b", "c", "d"))
  # a TOP loss recodes to itself, and the labels still come back
  rc2 <- frmtmb:::thres_pin_recode(p, list(resp_name = "y"), c(1, 2, 3),
                                   c("a", "b", "c"))
  expect_identical(rc2$y, c(1L, 2L, 3L))
  expect_identical(rc2$levels, c("a", "b", "c", "d"))
  # a category the fit never saw cannot be placed on its scale
  expect_error(
    frmtmb:::thres_pin_recode(p, list(resp_name = "y"), c(1, 2),
                              c("a", "z")),
    "never saw")
  # the ABSENT cases: no pin at all, and a response with no labels
  expect_identical(
    frmtmb:::thres_pin_recode(NULL, list(resp_name = "y"), c(1, 2), "a"),
    list(y = c(1, 2), levels = "a"))
  expect_identical(
    frmtmb:::thres_pin_recode(p, list(resp_name = "y"), c(1, 2), NULL),
    list(y = c(1, 2), levels = NULL))
})

test_that("the accepts_aterms half of the pin gate can be false", {
  skip_on_cran()
  # No family in the package makes it false, so the construction is an
  # ordinal family from outside that does not declare the term. Pinning
  # writes `thres`, which such a family would then refuse by name, so the
  # gate has to keep it out.
  shim <- function(at) {
    fam <- frmtmb_family("ordy", type = "ordinal", accepts_aterms = at,
                         dpars = "mu", links = list(mu = "identity"),
                         lpdf = function(y, dpars, aterms) 0 * y)
    list(frame = list(par_template = list(tau_raw = c(0, 0, 0)),
                      y_levels = list(y = c("1", "2", "3", "4"))),
         spec = list(responses = list(
           y = list(resp_name = "y", family = fam))))
  }
  expect_null(frmtmb:::thres_pin_of_fit(shim("weights")))
  p <- frmtmb:::thres_pin_of_fit(shim(c("weights", "thres")))
  expect_identical(p$y$nthres, 3L)
  expect_identical(p$y$levels_y, c("1", "2", "3", "4"))
})

test_that("a leave-one-out refit keeps an ordered factor's categories", {
  skip_on_cran()
  # the model frame drops a factor level that no row of the subset takes,
  # so the pinned count used to exceed the labels and the whole row of
  # the table came back NA
  dd <- refit_one_top()
  dd$y <- factor(dd$y, levels = 1:4, ordered = TRUE)
  itop <- which(as.integer(dd$y) == 4L)
  expect_length(itop, 1L)
  fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
  expect_length(fit$estimates[["tau_raw"]], 3L)
  inf <- suppressWarnings(influence(fit, force = TRUE))
  expect_false(anyNA(inf$fixed))
  expect_false(anyNA(cooks.distance(inf)))
  expect_identical(colnames(inf$fixed),
                   c("x", "tau_raw_1", "tau_raw_2", "tau_raw_3"))
  # the same row as the integer-coded fit gives, which is the point: the
  # spelling of the response must not change the diagnostic
  di <- refit_one_top()
  fi <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = di))
  ii <- suppressWarnings(influence(fi, force = TRUE))
  sp <- apply(ii$fixed, 2, stats::sd)
  expect_lt(max(abs(inf$fixed[itop, ] - ii$fixed[itop, ]) / sp), 1e-3)
})

test_that("a bootstrap replicate without the top category keeps the layout", {
  skip_on_cran()
  # the property frm_bootstrap() has because refit() reuses the
  # assembled frame rather than recounting; pinned so a future change to
  # refit() cannot lose it
  set.seed(202)
  x <- stats::rnorm(40)
  cp <- cbind(stats::plogis(-0.6 - 0.5 * x), stats::plogis(0.5 - 0.5 * x),
              stats::plogis(2.6 - 0.5 * x))
  dd <- data.frame(x = x, y = 1L + rowSums(stats::runif(40) > cp))
  FUNt <- function(f) {
    c(fixef(f, flatten = TRUE), tau = f$estimates[["tau_raw"]])
  }
  for (fam in list(cumulative(), sratio())) {
    fit <- suppressWarnings(frm(bf(y ~ x), family = fam, data = dd))
    set.seed(11)
    sims <- simulate(fit, nsim = 60, re_formula = NA)
    lost <- which(vapply(sims, function(v) max(as.integer(v)), 1L) < 4L)
    expect_gt(length(lost), 0L)
    bs <- suppressWarnings(frm_bootstrap(fit, FUN = FUNt, nsim = 60,
                                         seed = 11))
    expect_identical(dim(bs$t), c(60L, 4L))
    expect_identical(colnames(bs$t), c("x", "tau1", "tau2", "tau3"))
    expect_false(anyNA(bs$t))

    yb <- as.integer(sims[[lost[1L]]])
    rf <- suppressWarnings(refit(fit, yb))
    pinned <- suppressWarnings(
      frm(bf(y | thres(3) ~ x), family = fam,
          data = data.frame(x = dd$x, y = yb)))
    expect_length(rf$estimates[["tau_raw"]], 3L)
    ll1 <- as.numeric(logLik(rf))
    ll2 <- as.numeric(logLik(pinned))
    expect_lt(abs(ll1 - ll2) / abs(ll2), sqrt(.Machine$double.eps))
  }
})

test_that("a grouped bootstrap replicate keeps the per-level counts", {
  skip_on_cran()
  dd <- refit_grouped_top(503)
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                              family = cumulative(), data = dd))
  th <- fit$spec$responses$y$family[["thres"]]
  nraw <- length(fit$estimates[["tau_raw"]])
  set.seed(13)
  sims <- simulate(fit, nsim = 40, re_formula = NA)
  lost <- which(vapply(sims, function(v) {
    any(tapply(as.integer(v), dd$g, max) - 1L < th[["nthres"]])
  }, TRUE))
  expect_gt(length(lost), 0L)
  FUNt <- function(f) {
    c(fixef(f, flatten = TRUE), tau = f$estimates[["tau_raw"]])
  }
  bs <- suppressWarnings(frm_bootstrap(fit, FUN = FUNt, nsim = 40, seed = 13))
  expect_identical(dim(bs$t), c(40L, 1L + nraw))
  expect_false(anyNA(bs$t))

  yb <- as.integer(sims[[lost[1L]]])
  db <- data.frame(x = dd$x, g = dd$g, y = yb,
                   k = as.integer(th[["nthres"]][match(as.character(dd$g),
                                                       th[["levels"]])]))
  rf <- suppressWarnings(refit(fit, yb))
  pinned <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                                 family = cumulative(), data = db))
  expect_length(rf$estimates[["tau_raw"]], nraw)
  ll1 <- as.numeric(logLik(rf))
  ll2 <- as.numeric(logLik(pinned))
  expect_lt(abs(ll1 - ll2) / abs(ll2), sqrt(.Machine$double.eps))
})

test_that("a categorical refit keeps every category's predictor", {
  skip_on_cran()
  # the category count of categorical() is resolved on the SPEC, which
  # every in-package refit reuses, so the layout is safe there. The
  # complement of the ordinal case, checked rather than assumed.
  set.seed(504)
  x <- stats::rnorm(50)
  P <- cbind(1, exp(0.2 + 0.4 * x), exp(-2.6 + 0.3 * x))
  P <- P / rowSums(P)
  y <- apply(P, 1, function(p) sample.int(3L, 1L, prob = p))
  top <- which(y == 3L)
  y[top[-1L]] <- 2L
  dd <- data.frame(x = x, y = factor(y))
  expect_length(which(dd$y == "3"), 1L)
  fit <- suppressWarnings(frm(bf(y ~ x), family = categorical(), data = dd))
  inf <- suppressWarnings(influence(fit, force = TRUE))
  expect_identical(ncol(inf$fixed), 4L)
  expect_false(anyNA(inf$fixed))
  bs <- suppressWarnings(frm_bootstrap(fit, nsim = 10, seed = 21))
  expect_identical(dim(bs$t), c(10L, 4L))
  expect_false(anyNA(bs$t))
})

test_that("a prior draw never writes one value into a threshold vector", {
  skip_on_cran()
  # The recycling: one draw per ENTRY, and a class "Intercept" entry on
  # sratio, cratio and acat covers the whole unordered vector.
  e3 <- list(comp = "tau_raw", idx = 1:3, scale = "internal",
             dist = prior_normal(0, 2))
  expect_error(
    frmtmb:::draw_prior_pars(list(tau_raw = c(0, 0, 0)), list(e3),
                             "tau_raw[1:3]"),
    "density on 3 parameters at once")
  # the ABSENT case: one index is one draw for one parameter, so it must
  # still work
  e1 <- list(comp = "tau_raw", idx = 1L, scale = "internal",
             dist = prior_normal(0, 2))
  dr <- frmtmb:::draw_prior_pars(list(tau_raw = 0), list(e1), "tau_raw_1")
  expect_length(dr$est[["tau_raw"]], 1L)
  expect_true(is.finite(dr$est[["tau_raw"]]))

  # a label is one string, whatever the entry covers: the vector
  # spelling used to reach vapply() and die there
  set.seed(701)
  dd <- data.frame(x = stats::rnorm(40), y = rep(1:4, 10))
  bform <- frmtmb:::as_bform(bf(y ~ x), cratio())
  spec <- frmtmb:::parse_spec(bform)
  frame <- frmtmb:::assemble_frame(spec, dd)
  slots <- frmtmb:::nat_slots(frame, list(spec = spec, frame = frame,
                                          bform = bform))
  expect_length(frmtmb:::prior_entry_label(frame, slots, e3), 1L)
  expect_length(frmtmb:::prior_entry_label(frame, slots, e1), 1L)

  # and through the public path
  for (fam in list(cratio(), acat(), sratio())) {
    expect_error(
      frm_simulate(bf(y ~ x), family = fam, data = dd,
                   prior = set_prior("normal(0, 2)", class = "Intercept") +
                     set_prior("normal(0, 1)", class = "b"),
                   nsim = 2, seed = 702),
      "density on 3 parameters at once")
  }
  # cumulative keeps its own refusal, which names the ordering
  expect_error(
    frm_simulate(bf(y ~ x), family = cumulative(), data = dd,
                 prior = set_prior("normal(0, 2)", class = "Intercept") +
                   set_prior("normal(0, 1)", class = "b"),
                 nsim = 2, seed = 702),
    "would not produce an ordered one")
  # two categories is one threshold, and the draw must still happen
  set.seed(703)
  dd2 <- data.frame(x = stats::rnorm(40), y = rep(1:2, 20))
  s <- frm_simulate(bf(y ~ x), family = cratio(), data = dd2,
                    prior = set_prior("normal(0, 2)", class = "Intercept") +
                      set_prior("normal(0, 1)", class = "b"),
                    nsim = 2, seed = 704)
  p <- attr(s, "pars")
  expect_identical(names(p), c("tau_raw_1", "b_x"))
  expect_true(all(is.finite(unlist(p))))
})

# The three response codings an ordinal fit accepts. An unordered factor
# is refused by frm() itself, in brms's words, so it is not here.
refit_codings <- list(
  integer = function(v) v,
  ordered = function(v) factor(v, levels = 1:4, ordered = TRUE),
  character = function(v) as.character(v))

# n = 60, four categories, EXACTLY ONE row in the INTERIOR category 2,
# so deleting it empties a category that is not the top one.
refit_one_interior <- function(seed = 901, n = 60) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(0.9 - 0.5 * x), stats::plogis(1.0 - 0.5 * x),
              stats::plogis(2.0 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  i2 <- which(y == 2L)
  if (!length(i2)) y[1L] <- 2L else y[i2[-1L]] <- 1L
  data.frame(x = x, y = y)
}

test_that("an emptied INTERIOR category does not slide the row over", {
  skip_on_cran()
  # `drop.unused.levels = TRUE` removes the emptied level from the model
  # frame, which renumbers every category above it, so a factor response
  # used to report tau_raw_3's value in tau_raw_2's column with nothing
  # but a trailing NA to show it.
  dd <- refit_one_interior()
  i2 <- which(dd$y == 2L)
  expect_length(i2, 1L)
  got <- list()
  for (nm in names(refit_codings)) {
    d <- data.frame(x = dd$x, y = refit_codings[[nm]](dd$y))
    fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d))
    expect_length(fit$estimates[["tau_raw"]], 3L)
    inf <- suppressWarnings(influence(fit, force = TRUE))
    expect_false(anyNA(inf$fixed))
    expect_false(anyNA(cooks.distance(inf)))
    got[[nm]] <- inf$fixed
  }
  # bitwise, on the whole table: the coding of the response must not
  # change the diagnostic at all
  expect_identical(got$ordered, got$integer)
  expect_identical(got$character, got$integer)

  # and the emptied category's threshold is the unidentified one, in its
  # OWN column, matching a hand-pinned fit everywhere else
  ref <- suppressWarnings(
    frm(bf(y | thres(3) ~ x), family = cumulative(),
        data = dd[-i2, , drop = FALSE]))
  want <- frmtmb:::get_coef.frmtmb_fit(ref)
  row <- got$integer[i2, ]
  expect_identical(names(want), colnames(got$integer))
  ident <- setdiff(names(want), "tau_raw_2")
  spread <- apply(got$integer, 2, stats::sd)
  expect_lt(max(abs(row[ident] - want[ident]) / spread[ident]), 1e-3)
  scale_id <- max(abs(want[c("tau_raw_1", "tau_raw_3")]))
  expect_gt(abs(row[["tau_raw_2"]]) / scale_id, 3)
  expect_gt(abs(want[["tau_raw_2"]]) / scale_id, 3)
})

test_that("an emptied category does not slide a grouped row over", {
  skip_on_cran()
  # the level only disappears when the category is absent from the WHOLE
  # response, not merely from one thres(gr = ) level
  set.seed(904)
  n <- 96
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  y[y == 2L] <- 1L
  y[which(g == "a")[1L]] <- 2L
  dd <- data.frame(x = x, g = g, y = y)
  ih <- which(dd$y == 2L)
  expect_length(ih, 1L)

  got <- list()
  for (nm in names(refit_codings)) {
    d <- data.frame(x = dd$x, g = dd$g, y = refit_codings[[nm]](dd$y))
    fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                                family = cumulative(), data = d))
    inf <- suppressWarnings(influence(fit, force = TRUE))
    expect_false(anyNA(inf$fixed))
    got[[nm]] <- inf$fixed
  }
  expect_identical(got$ordered, got$integer)
  expect_identical(got$character, got$integer)
  # the shift used to run over five coefficients with three trailing NA,
  # so what the assertion is about is POSITION
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                              family = cumulative(), data = dd))
  th <- fit$spec$responses$y$family[["thres"]]
  ds <- dd[-ih, , drop = FALSE]
  ds$k <- as.integer(th[["nthres"]][match(as.character(ds$g),
                                          th[["levels"]])])
  ref <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                              family = cumulative(), data = ds))
  want <- frmtmb:::get_coef.frmtmb_fit(ref)
  row <- got$integer[ih, ]
  expect_identical(names(want), colnames(got$integer))
  # the unidentified thresholds run off toward infinity on both sides;
  # every other coefficient has to line up slot for slot
  spread <- apply(got$integer, 2, stats::sd)
  big <- names(want)[abs(want) > 10]
  ident <- setdiff(names(want), big)
  expect_lt(max(abs(row[ident] - want[ident]) / spread[ident]), 1e-3)
  # level a's second threshold bounds the emptied category
  expect_true("tau_raw_2" %in% big)
  expect_gt(length(ident), 3L)
})

test_that("a hand-written thres(K) is refit on every response coding", {
  skip_on_cran()
  # a user-written count used to short-circuit the recoding as well as
  # the pinned count, so the ordered factor was the one coding whose
  # deleted row came back entirely NA
  dd <- refit_one_top()
  itop <- which(dd$y == 4L)
  got <- list()
  for (nm in names(refit_codings)) {
    d <- data.frame(x = dd$x, y = refit_codings[[nm]](dd$y))
    fit <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = cumulative(),
                                data = d))
    inf <- suppressWarnings(influence(fit, force = TRUE))
    expect_false(anyNA(inf$fixed))
    got[[nm]] <- inf$fixed
  }
  expect_identical(got$ordered, got$integer)
  expect_identical(got$character, got$integer)
  expect_false(anyNA(got$ordered[itop, ]))
})

test_that("influence() counts the refits that failed and says so", {
  skip_on_cran()
  # A table of nothing but NA used to come back without a word. Every
  # construction here is a refit that CANNOT represent the fitted model.
  set.seed(1101)
  n <- 50
  x <- stats::rnorm(n)
  y <- 1L + rowSums(cbind(stats::runif(n) > stats::plogis(-0.4 - 0.5 * x),
                          stats::runif(n) > stats::plogis(1.2 - 0.5 * x)))
  dd <- data.frame(x = x, y = pmin(y, 3L))
  fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
  expect_length(fit$estimates[["tau_raw"]], 2L)

  # every unit fails: two rows of `data` reach a category the fit never
  # saw, so no single deletion removes both
  d2 <- dd
  d2$y[1:2] <- 4L
  expect_error(suppressWarnings(influence(fit, data = d2, force = TRUE)),
               "all 50 deletion refits failed")

  # exactly one row does, so one unit refits and the rest are NA: that is
  # a warning carrying the count, and the table still comes back
  d3 <- dd
  d3$y[7L] <- 4L
  expect_warning(inf <- influence(fit, data = d3, force = TRUE),
                 "49 of 50 deletion refits failed")
  expect_identical(unname(which(apply(inf$fixed, 1,
                                      function(z) !anyNA(z)))), 7L)

  # the ABSENT case: nothing wrong warns about nothing
  expect_silent(suppressMessages(influence(fit, force = TRUE)))
})

test_that("a refit cannot add or drop a thres(gr = ) level", {
  skip_on_cran()
  set.seed(1103)
  n <- 96
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  dd <- data.frame(x = x, g = g, y = y)
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                              family = cumulative(), data = dd))

  # a level the fit never saw
  d2 <- dd
  levels(d2$g) <- c(levels(dd$g), "d")
  d2$g[1:4] <- "d"
  expect_error(suppressWarnings(influence(fit, data = d2, force = TRUE)),
               "cannot add a level")

  # a groups = deletion that empties a whole level: its thresholds have
  # no data at all, so the unit is refused rather than dropped, which
  # used to slide every later level's thresholds one column early
  d3 <- dd
  d3$idg <- d3$g
  fit3 <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | idg)),
                               family = cumulative(), data = d3))
  expect_error(suppressWarnings(influence(fit3, groups = "idg")),
               "cannot drop a level")

  # the ABSENT case, and the common one: a grouping factor that is NOT
  # the threshold grouping factor leaves every level populated
  set.seed(1201)
  n4 <- 120
  g4 <- factor(rep(c("a", "b", "c"), length.out = n4))
  id4 <- factor(rep(1:8, length.out = n4))
  x4 <- stats::rnorm(n4)
  y4 <- vapply(seq_len(n4), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g4[i])]] - 0.5 * x4[i]))
  }, 1L)
  d4 <- data.frame(x = x4, g = g4, id = id4, y = y4)
  fit4 <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | id)),
                               family = cumulative(), data = d4))
  inf4 <- influence(fit4, groups = "id")
  expect_false(anyNA(inf4$fixed))
  expect_false(anyNA(cooks.distance(inf4)))
})
