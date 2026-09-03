# Several condition sets are FACETS of one page, not one page each, and
# each facet carries only its own condition's raw observations. The
# display is checked through the structures it computes (the panel grid
# and the per-condition point subsets), not through pixels.

skip_on_cran()

sim_facet <- function(n = 120, seed = 41) {
  set.seed(seed)
  dd <- data.frame(
    x = stats::rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(seq_len(12), each = n / 12))
  )
  dd$y <- 1 + 0.8 * dd$x + as.integer(dd$f) + stats::rnorm(n, 0, 0.5)
  dd
}

facet_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dd <- sim_facet()
      cache <<- list(dd = dd,
                     fit = frm(bf(y ~ x + f + (1 | g)),
                               family = gaussian(), data = dd))
    }
    cache
  }
})

f_conditions <- function(dd) {
  data.frame(f = factor(c("a", "b", "c"), levels = levels(dd$f)),
             row.names = c("f = a", "f = b", "f = c"))
}

test_that("the panel grid is roughly square and honors ncol", {
  # brms passes ncol straight to facet_wrap(), whose NULL default is
  # roughly square
  expect_identical(ce_facet_layout(1L), c(nrow = 1L, ncol = 1L))
  expect_identical(ce_facet_layout(4L), c(nrow = 2L, ncol = 2L))
  expect_identical(ce_facet_layout(9L), c(nrow = 3L, ncol = 3L))
  expect_identical(ce_facet_layout(10L), c(nrow = 3L, ncol = 4L))
  expect_identical(ce_facet_layout(10L, 5L), c(nrow = 2L, ncol = 5L))
  expect_identical(ce_facet_layout(10L, 1L), c(nrow = 10L, ncol = 1L))
  # every panel always fits on the one page
  for (n in 1:12) {
    for (nc in list(NULL, 1L, 3L, 5L)) {
      lay <- ce_facet_layout(n, nc)
      expect_gte(lay[["nrow"]] * lay[["ncol"]], n)
    }
  }
  # asking for more columns than there are panels only wastes the page
  expect_identical(ce_facet_layout(3L, 9L), c(nrow = 1L, ncol = 3L))
  expect_error(plot(structure(list(), class = "frmtmb_conditional_effects"),
                    ncol = 0), "ncol")
})

test_that("each condition's panel gets only its own observations", {
  cs <- facet_fit()
  ce <- conditional_effects(cs$fit, effects = "x", resolution = 8,
                            conditions = f_conditions(cs$dd))
  pts <- attr(ce[["x"]], "points_df")
  expect_true("cond__" %in% names(pts))
  expect_identical(levels(pts$cond__), c("f = a", "f = b", "f = c"))
  # one condition per level of a factor: every observation is drawn
  # exactly once, in its own panel
  expect_identical(nrow(pts), nrow(cs$dd))
  for (lv in levels(pts$cond__)) {
    got <- pts[pts$cond__ == lv, ]
    ref <- cs$dd[cs$dd$f == sub("f = ", "", lv), ]
    expect_identical(sort(got$y), sort(ref$y))
    expect_identical(sort(got$x), sort(ref$x))
  }
  # a single condition set is not faceted, so its points are not split
  ce1 <- conditional_effects(cs$fit, effects = "x", resolution = 8)
  expect_null(attr(ce1[["x"]], "points_df")[["cond__"]])
})

test_that("a numeric condition keeps every observation in every panel", {
  # brms's select_points = 0 rule, measured against its make_point_frame:
  # a reference value on a continuum names no observation, so filtering
  # on it would empty the panels rather than subset them
  cs <- facet_fit()
  cn <- data.frame(x = c(-1, 0, 1), row.names = c("x=-1", "x=0", "x=1"))
  ce <- conditional_effects(cs$fit, effects = "f", resolution = 8,
                            conditions = cn)
  pts <- attr(ce[["f"]], "points_df")
  expect_identical(as.integer(table(pts$cond__)),
                   rep(nrow(cs$dd), 3L))

  # a grouping factor is a label, not a continuum, so it still matches
  cg <- data.frame(g = factor(c("1", "2"), levels = levels(cs$dd$g)),
                   row.names = c("g1", "g2"))
  ceg <- conditional_effects(cs$fit, effects = "x", resolution = 8,
                             conditions = cg)
  ptg <- attr(ceg[["x"]], "points_df")
  expect_identical(as.integer(table(ptg$cond__)),
                   as.integer(table(cs$dd$g))[1:2])
})

test_that("the numeric layer is untouched by the display change", {
  cs <- facet_fit()
  cnd <- f_conditions(cs$dd)
  ce <- conditional_effects(cs$fit, effects = "x", resolution = 8,
                            conditions = cnd)
  df <- ce[["x"]]
  # the effect data frame gains no column for the display: the panels
  # come from cond__, which was already there
  expect_identical(names(df),
                   c("x", "estimate__", "se__", "lower__", "upper__",
                     "cond__"))
  expect_identical(nrow(df), 8L * nrow(cnd))
  # the curves are the per-condition curves, unchanged by faceting
  for (i in seq_len(nrow(cnd))) {
    one <- conditional_effects(cs$fit, effects = "x", resolution = 8,
                               conditions = as.list(cnd[i, , drop = FALSE]))
    sub <- df[df$cond__ == rownames(cnd)[i], ]
    expect_equal(sub$estimate__, one[["x"]]$estimate__)
    expect_equal(sub$lower__, one[["x"]]$lower__)
    expect_equal(sub$upper__, one[["x"]]$upper__)
  }
  # drawing does not mutate the object it was handed
  before <- ce
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(ce, points = TRUE, ask = FALSE)
  expect_identical(ce, before)
})

test_that("both faceted paths draw one page and leave par alone", {
  cs <- facet_fit()
  ce <- conditional_effects(cs$fit, effects = c("x", "f"), resolution = 8,
                            conditions = f_conditions(cs$dd))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  mfrow0 <- graphics::par("mfrow")
  rng <- .Random.seed
  expect_no_error(plot(ce, points = TRUE, ncol = 2, ask = FALSE))
  # the deterministic point spread must not disturb the user's stream
  expect_identical(.Random.seed, rng)
  # a faceted page that left mfrow set would swallow the next plot
  expect_identical(graphics::par("mfrow"), mfrow0)
  expect_no_error(plot(ce, ask = FALSE))
  expect_no_error(plot(ce, points = TRUE, ncol = 1, ask = FALSE))

  # the fallback grid: same page, same ncol, no tinyplot
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "tinyplot")) FALSE
      else base::requireNamespace(package, ...)
    }, .package = "base")
  expect_no_error(plot(ce, points = TRUE, ncol = 2, ask = FALSE))
  expect_identical(graphics::par("mfrow"), mfrow0)
})

test_that("a two-predictor effect and an ordinal display both facet", {
  cs <- facet_fit()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # the second predictor takes the grouping slot, cond__ the facet slot
  cn <- data.frame(x = c(-1, 1), row.names = c("lo", "hi"))
  ce2 <- conditional_effects(cs$fit, effects = "x:f", resolution = 8,
                             conditions = cn)
  expect_true("cond__" %in% names(ce2[["x:f"]]))
  expect_no_error(plot(ce2, points = TRUE, ncol = 2, ask = FALSE))

  # the per-category ordinal display already spends both slots, so it
  # takes the base grid rather than a second facet dimension
  set.seed(7)
  od <- data.frame(x = stats::rnorm(90),
                   f = factor(rep(c("a", "b"), length.out = 90)))
  od$y <- factor(cut(od$x + stats::rnorm(90), 3), labels = 1:3,
                 ordered = TRUE)
  ofit <- frm(bf(y ~ x + f), family = cumulative(), data = od)
  oce <- conditional_effects(
    ofit, effects = "x", resolution = 6,
    conditions = data.frame(f = factor(c("a", "b"), levels = c("a", "b")),
                            row.names = c("f = a", "f = b")))
  expect_true("cats__" %in% names(oce[["x"]]))
  expect_no_error(plot(oce, ncol = 2, ask = FALSE))
})
