# plot() of a conditional-effects object and of a hypothesis takes the
# arguments of brms 2.23.0's plot.brms_conditional_effects() and
# plot.brmshypothesis(), drawing with base graphics. On 0.66.0 every one
# of these calls stopped at "plot() has no argument" (the ported rows
# brmsfit-methods:154, :162, :164, :391 and :396;
# dev/ceplot-log/probe-base.txt).

ce_plot_data <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(5)
      dd <- data.frame(x = stats::rnorm(80), z = stats::rnorm(80),
                       f = factor(rep(c("a", "b"), 40)))
      dd$y <- stats::rnorm(80, 1 + 0.5 * dd$x + (dd$f == "b") +
                             0.3 * dd$x * dd$z, 1)
      fit <- frm(bf(y ~ x * f + x * z), family = gaussian(), data = dd)
      cache <<- list(fit = fit, ce = conditional_effects(
        fit, effects = c("x", "f"), resolution = 5))
    }
    cache
  }
})

# how many pages a call draws: base graphics runs the "before.plot.new"
# hook for every new frame, and par("page") says whether that frame
# starts a new page (a panel of a par(mfrow) grid does not)
pages <- function(expr) {
  n <- 0L
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  old <- getHook("before.plot.new")
  on.exit(setHook("before.plot.new", old, "replace"), add = TRUE)
  setHook("before.plot.new", function() {
    if (isTRUE(graphics::par("page"))) n <<- n + 1L
  })
  force(expr)
  n
}

# the graphics operations one page records, by C entry point: a line is
# C_plotXY (as is the empty frame plot() opens), a rug is C_axis
display_ops <- function(expr) {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  grDevices::dev.control("enable")
  force(expr)
  rp <- grDevices::recordPlot()
  vapply(rp[[1L]], function(e) {
    f <- e[[2L]][[1L]]
    if (is.list(f) && !is.null(f$name)) f$name else ""
  }, "")
}

test_that("plot() takes every argument brms's two plot methods take", {
  # brms 2.23.0's formals, from dev/ceplot-brmssrc/; kept here rather
  # than read from brms so that a machine without brms checks them too
  brms_ce <- c("x", "ncol", "points", "rug", "mean", "jitter_width",
               "stype", "line_args", "cat_args", "errorbar_args",
               "surface_args", "spaghetti_args", "point_args",
               "rug_args", "facet_args", "theme", "ask", "plot", "...")
  brms_hyp <- c("x", "nvariables", "N", "ignore_prior", "chars",
                "colors", "theme", "ask", "plot", "...")
  ce_f <- names(formals(getS3method("plot",
                                    "frmtmb_conditional_effects")))
  hyp_f <- names(formals(getS3method("plot", "frmtmb_hypothesis")))
  expect_identical(ce_f, brms_ce)
  expect_identical(hyp_f, brms_hyp)
  skip_if_not_installed("brms")
  expect_identical(
    names(formals(utils::getFromNamespace("plot.brms_conditional_effects",
                                          "brms"))), brms_ce)
  expect_identical(
    names(formals(utils::getFromNamespace("plot.brmshypothesis", "brms"))),
    brms_hyp)
})

test_that("plot = FALSE returns one plot object per effect, undrawn", {
  s <- ce_plot_data()
  expect_identical(pages(p <- plot(s$ce, plot = FALSE)), 0L)
  expect_named(p, names(s$ce))
  expect_true(all(vapply(p, inherits, NA, "frmtmb_ce_plot")))
  # printing an object draws it, one page, as printing a ggplot does
  expect_identical(pages(print(p$x)), 1L)
  expect_identical(pages(plot(p$f)), 1L)
  # the list comes back invisibly whether or not the call draws
  expect_false(withVisible(plot(s$ce, plot = FALSE))$visible)
  expect_identical(pages(v <- withVisible(plot(s$ce, ask = FALSE))), 2L)
  expect_false(v$visible)
  expect_named(v$value, names(s$ce))
  # brms's deprecated alias, with brms's warning
  n <- allow_warnings(pages(plot(s$ce, do_plot = FALSE)),
                      "'do_plot' is deprecated",
                      require = "'do_plot' is deprecated")
  expect_identical(n, 0L)
})

test_that("rug draws the observed values of a numeric predictor", {
  s <- ce_plot_data()
  x_plain <- sum(display_ops(print(plot(s$ce, plot = FALSE)$x)) == "C_axis")
  x_rug <- sum(display_ops(print(plot(s$ce, plot = FALSE,
                                      rug = TRUE)$x)) == "C_axis")
  expect_identical(x_rug - x_plain, 1L)
  # a top rug too; a left or right one has nothing to mark, as in brms
  x_bt <- sum(display_ops(print(plot(s$ce, plot = FALSE, rug = TRUE,
                                     rug_args = list(sides = "btl"))$x)) ==
                "C_axis")
  expect_identical(x_bt - x_plain, 2L)
  # a factor predictor has no rug
  f_plain <- sum(display_ops(print(plot(s$ce, plot = FALSE)$f)) == "C_axis")
  f_rug <- sum(display_ops(print(plot(s$ce, plot = FALSE,
                                      rug = TRUE)$f)) == "C_axis")
  expect_identical(f_rug, f_plain)
})

test_that("mean = FALSE leaves the estimate out only beside spaghetti", {
  s <- ce_plot_data()
  lines_of <- function(ce, ...) {
    sum(display_ops(print(plot(ce, plot = FALSE, ...)$x)) == "C_plotXY")
  }
  # without spaghetti the estimate line IS the display, as in brms
  expect_identical(lines_of(s$ce, mean = FALSE), lines_of(s$ce))
  sp <- conditional_effects(s$fit, "x", resolution = 5, band = "boot",
                            boot = 5, seed = 1, spaghetti = TRUE)
  expect_identical(lines_of(sp) - lines_of(sp, mean = FALSE), 1L)
})

test_that("stype picks contour lines or a raster for a surface", {
  s <- ce_plot_data()
  su <- conditional_effects(s$fit, "x:z", surface = TRUE, resolution = 6)
  ops_c <- display_ops(print(plot(su, plot = FALSE)[[1L]]))
  ops_r <- display_ops(print(plot(su, plot = FALSE,
                                  stype = "raster")[[1L]]))
  expect_true("C_contour" %in% ops_c)
  expect_false("C_image" %in% ops_c)
  expect_true("C_image" %in% ops_r)
  expect_false("C_contour" %in% ops_r)
  expect_error(plot(su, stype = "tile"), "stype")
})

test_that("layer arguments translate, and the rest are named in a warning", {
  s <- ce_plot_data()
  expect_no_warning(pages(plot(
    s$ce, ask = FALSE, points = TRUE, rug = TRUE,
    line_args = list(colour = "red", linewidth = 2, fill = "orange",
                     alpha = 0.3),
    cat_args = list(size = 2, shape = 17),
    errorbar_args = list(width = 0.3, colour = "darkgreen"),
    point_args = list(size = 1, width = 0.1),
    rug_args = list(sides = "b", colour = "purple"),
    spaghetti_args = list(colour = "blue"),
    surface_args = list(bins = 5),
    facet_args = list(ncol = 1))))
  n <- allow_warnings(
    pages(plot(s$ce, ask = FALSE, line_args = list(stat = "identity"),
               cat_args = list(position = "dodge"))),
    "plot() ignores",
    require = c("line_args$stat", "cat_args$position"))
  expect_identical(n, 2L)
  # brms's own refusals of a malformed list
  expect_error(plot(s$ce, line_args = list(1)),
               "Argument 'line_args' must be named.", fixed = TRUE)
  expect_error(plot(s$ce, point_args = list(mapping = 1)),
               "Argument(s) mapping cannot be replaced.", fixed = TRUE)
  expect_error(plot(s$ce, rug_args = list(sides = "x")), "sides")
  expect_error(plot(s$ce, facet_args = list(scales = "loose")),
               "facet_args$scales", fixed = TRUE)
  expect_error(plot(s$ce, foo = 1), "no argument `foo`")
})

test_that("theme, jitter_width and points follow brms's rules", {
  s <- ce_plot_data()
  expect_error(plot(s$ce, theme = "bw"),
               "Argument 'theme' should be a 'theme' object.", fixed = TRUE)
  n <- allow_warnings(pages(plot(s$ce, ask = FALSE, jitter_width = 0.1,
                                 points = TRUE)),
                      "'jitter_width' is deprecated",
                      require = "'jitter_width' is deprecated")
  expect_identical(n, 2L)
  set.seed(24)
  d <- data.frame(x = stats::runif(80))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(80, 0, 0.3)
  cs <- conditional_smooths(frm(bf(y ~ s(x)), family = gaussian(),
                                data = d))
  expect_error(plot(cs, points = TRUE),
               paste0("Argument 'points' is invalid for objects returned ",
                      "by 'conditional_smooths'."), fixed = TRUE)
  expect_identical(pages(plot(cs, ask = FALSE)), 1L)
  skip_if_not_installed("ggplot2")
  n <- allow_warnings(pages(plot(s$ce, ask = FALSE,
                                 theme = ggplot2::theme_bw())),
                      "is a ggplot2 theme",
                      require = "is a ggplot2 theme")
  expect_identical(n, 2L)
})

test_that("facet_args lays out the panels of several conditions", {
  expect_identical(ce_facet_layout(4L, nrow = 1L), c(nrow = 1L, ncol = 4L))
  expect_identical(ce_facet_layout(4L, ncol = 1L, nrow = 3L),
                   c(nrow = 4L, ncol = 1L))
  s <- ce_plot_data()
  cc <- conditional_effects(s$fit, "x", resolution = 5,
                            conditions = data.frame(f = c("a", "b")))
  expect_identical(pages(plot(cc, ask = FALSE,
                              facet_args = list(nrow = 1L,
                                                scales = "free"))), 1L)
})

test_that("plot() of a hypothesis takes brms's arguments", {
  s <- ce_plot_data()
  h <- hypothesis(s$fit, c("x > 0", "fb = 0", "x + z = 0"))
  # nvariables hypotheses to a page, stacked, as brms's facet_wrap()
  expect_identical(pages(plot(h, ask = FALSE)), 1L)
  expect_identical(pages(p <- plot(h, nvariables = 2, plot = FALSE)), 0L)
  expect_length(p, 2L)
  expect_true(all(vapply(p, inherits, NA, "frmtmb_hyp_plot")))
  expect_identical(pages(print(p[[2L]])), 1L)
  expect_identical(pages(plot(h, nvariables = 2, ask = FALSE)), 2L)
  n <- allow_warnings(pages(plot(h, N = 1, ask = FALSE)),
                      "Argument 'N' is deprecated",
                      require = "Argument 'N' is deprecated")
  expect_identical(n, 3L)
  # a fit has no prior draws, so there is nothing for ignore_prior to
  # leave out; brms draws the posterior alone then too
  expect_identical(pages(plot(h, ignore_prior = TRUE, ask = FALSE)), 1L)
  expect_error(plot(h, colors = c("red", "blue", "green")),
               "Argument 'colors' must be of length 2.", fixed = TRUE)
  expect_error(plot(h, theme = 1),
               "Argument 'theme' should be a 'theme' object.", fixed = TRUE)
  expect_false(withVisible(plot(h, plot = FALSE))$visible)
})

test_that("a long hypothesis is cut in its title as brms cuts it", {
  # brms 2.23.0's limit_chars() on the same labels
  # (dev/ceplot-limitchars.R)
  x <- c("(x+fb+z+x:fb+x:z)-(0.123456789) = 0", "(x) > 0",
         "(abcdefghijklmnopqrst) = 0", "(abcdefghijklmnopq) < 0")
  expect_identical(hyp_limit_chars(x, 20),
                   c("(x+fb+z+x:fb+x:z)... = 0", "(x) > 0",
                     "(abcdefghijklmnop... = 0", "(abcdefghijklmnopq) < 0"))
  expect_identical(hyp_limit_chars(x, 5),
                   c("(x... = 0", "(x) > 0", "(a... = 0", "(a... < 0"))
  expect_identical(hyp_limit_chars(x, NULL), x)
})

test_that("graphical parameters are named when they are ignored", {
  # review m7: col, main, xlab and lwd passed in silence, against the
  # rule that a setting the caller asked for must not
  s <- ce_plot_data()
  n <- allow_warnings(pages(plot(s$ce, ask = FALSE, col = "red",
                                 main = "m")),
                      "plot() ignores the graphical parameter(s)",
                      require = "`col`, `main`")
  expect_identical(n, 2L)
  h <- hypothesis(s$fit, "x > 0")
  n <- allow_warnings(pages(plot(h, ask = FALSE, xlab = "v", lwd = 2)),
                      "plot() ignores the graphical parameter(s)",
                      require = "`xlab`, `lwd`")
  expect_identical(n, 1L)
  p <- plot(s$ce, plot = FALSE)
  allow_warnings(pages(plot(p$x, col = 2)), "plot() ignores",
                 require = "`col`")
  # the refusal of an unknown argument does not list the hidden alias
  e <- tryCatch(plot(s$ce, foo = 1), error = function(e) conditionMessage(e))
  expect_match(e, "no argument `foo`", fixed = TRUE)
  expect_false(grepl("do_plot", e, fixed = TRUE))
  # ask = NULL, the default of 0.66.0, is taken as TRUE
  expect_identical(pages(plot(s$ce, ask = NULL)), 2L)
  expect_identical(pages(plot(h, ask = NULL)), 1L)
})

test_that("printing a plot object treats graphical parameters as plot()", {
  # punch round 2: print(p$x, main = ) stopped where plot(p$x, main = )
  # warned
  s <- ce_plot_data()
  p <- plot(s$ce, plot = FALSE)
  n <- allow_warnings(pages(print(p$x, main = "t")), "plot() ignores",
                      require = "`main`")
  expect_identical(n, 1L)
  h <- hypothesis(s$fit, "x > 0")
  ph <- plot(h, plot = FALSE)
  n <- allow_warnings(pages(print(ph[[1L]], main = "t")), "plot() ignores",
                      require = "`main`")
  expect_identical(n, 1L)
  # `col` abbreviates `colors` and is that argument, as in brms
  expect_error(plot(h, col = "red"),
               "Argument 'colors' must be of length 2.", fixed = TRUE)
})
