# plot.frmtmb_influence: the Cook's distance index plot and the dfbetas
# panels, on both deletion kinds. The fixtures are the ones the existing
# influence tests use, so the planted outlier is a known answer.

# the group-deletion fixture of test-v15.R: group 7 gets an extra slope
sim_infl_groups <- function(seed = 23) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
  dd$y <- 1 + 0.5 * dd$x + rnorm(12, 0, 0.4)[dd$g] + rnorm(120, 0, 0.5)
  dd$y[dd$g == "7"] <- dd$y[dd$g == "7"] + 4 * dd$x[dd$g == "7"]
  dd
}

test_that("plot() labels the planted outlier group and bands the dfbetas", {
  dd <- sim_infl_groups()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  infl <- influence(fit, groups = "g")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  drawn <- character(0)
  hlines <- numeric(0)
  testthat::local_mocked_bindings(
    text = function(x, y, labels, ...) {
      drawn <<- c(drawn, as.character(labels))
      invisible(NULL)
    },
    abline = function(h = NULL, ...) {
      hlines <<- c(hlines, if (is.null(h)) numeric(0) else h)
      invisible(NULL)
    },
    .package = "graphics"
  )

  # panel 1 alone: the most influential case is the planted group
  expect_identical(plot(infl, which = 1, ask = FALSE, labels = 1), infl)
  expect_identical(drawn, "7")

  # every panel: Cook's distance plus one dfbetas panel per coefficient
  drawn <- character(0)
  hlines <- numeric(0)
  expect_no_error(plot(infl, ask = FALSE, labels = 1))
  np <- 1L + ncol(dfbetas(infl))
  expect_length(drawn, np)
  expect_identical(drawn[1L], "7")

  # the conventional +/- 2 / sqrt(n) cutoff, on the number of deleted
  # units, is drawn on every dfbetas panel and on no other
  band <- 2 / sqrt(nrow(infl$fixed))
  expect_equal(sort(unique(hlines)), c(-band, 0, band))
  expect_length(hlines, np + 2L * (np - 1L))

  # more labels than the default, and the top case stays first
  drawn <- character(0)
  plot(infl, which = 1, ask = FALSE)
  expect_length(drawn, 3L)
  expect_identical(drawn[1L], "7")
})

test_that("which= selects panels and an all-NA table is refused", {
  dd <- sim_infl_groups()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  infl <- influence(fit, groups = "g")
  cn <- colnames(dfbetas(infl))

  panels <- list()
  local_mocked_bindings(
    infl_index_panel = function(v, units, xlab, ylab, main, labels,
                               band = NULL, ...) {
      panels[[length(panels) + 1L]] <<- list(main = main, xlab = xlab,
                                             band = band)
      invisible(NULL)
    }
  )
  # 1 + j is the jth coefficient's dfbetas panel
  plot(infl, which = 2, ask = FALSE)
  expect_identical(vapply(panels, `[[`, "", "main"), cn[1L])
  expect_equal(panels[[1L]]$band, 2 / sqrt(nrow(infl$fixed)))

  panels <- list()
  plot(infl, which = c(1L, 3L), ask = FALSE)
  expect_identical(vapply(panels, `[[`, "", "main"),
                   c("Cook's distance", cn[2L]))
  # the Cook's panel carries no reference band: there is no conventional
  # cutoff for it, unlike dfbetas
  expect_null(panels[[1L]]$band)
  expect_true(all(vapply(panels, `[[`, "", "xlab") == "Level of 'g'"))

  panels <- list()
  expect_no_error(plot(infl, which = integer(0), ask = FALSE))
  expect_length(panels, 0L)
})

test_that("a table of failed refits is refused rather than plotted", {
  dd <- sim_infl_groups()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  infl <- influence(fit, groups = "g")
  infl$fixed[] <- NA_real_
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(infl, ask = FALSE), "Every deletion refit failed")
})

test_that("plot() runs on an observation-deletion object", {
  set.seed(11)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 0.5)
  dd$y[5] <- dd$y[5] + 6
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  infl <- influence(fit)
  expect_null(infl$groups)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  drawn <- character(0)
  testthat::local_mocked_bindings(
    text = function(x, y, labels, ...) {
      drawn <<- c(drawn, as.character(labels))
      invisible(NULL)
    },
    .package = "graphics"
  )
  expect_no_error(plot(infl, ask = FALSE, labels = 1))
  expect_identical(drawn[1L], "5")

  # the deleted unit is a row here, and the axis label says so
  xlabs <- character(0)
  local_mocked_bindings(
    infl_index_panel = function(v, units, xlab, ...) {
      xlabs <<- c(xlabs, xlab)
      invisible(NULL)
    }
  )
  plot(infl, ask = FALSE)
  expect_true(all(xlabs == "Observation"))
})
