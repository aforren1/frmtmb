# pp_check() on a fit, for every type bayesplot exposes, against what
# brms 2.23.0 does on the same designs (dev/correct-log/brms-ppcheck.txt
# and brms-ppcheck2.txt, script dev/correct-brms-probe.R). Through 0.61.0
# the fit method passed `group = "g"` to bayesplot as one string, so
# every *_grouped type died on "length(group) must be equal to the
# number of observations", `x` the same way on "is.numeric(x) is not
# TRUE", and type = "violin" on "object 'ppc_violin' not found".

skip_if_not_installed("bayesplot")
skip_if_not_installed("ggplot2")

ppt_data <- function() {
  set.seed(110)
  n <- 60
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(letters[1:6], length.out = n)))
  eta <- 0.3 + 0.5 * d$x + stats::rnorm(6, 0, 0.5)[d$g]
  d$y <- eta + stats::rnorm(n)
  d$cnt <- stats::rpois(n, exp(eta))
  d$bin <- stats::rbinom(n, 1, stats::plogis(eta - 0.3))
  d
}

# One gaussian fit for the many one-behavior blocks below. Each refusal
# gets its own test_that, because testthat ENDS a block at an error and
# 0.61.0 errors on the first of them, which would hide the rest on the
# base arm; refitting per block would be the only cost of that, so the
# fit is made once and handed out.
ppt_gauss_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- frm(bf(y ~ x + (1 | g)) + gaussian(), data = ppt_data())
    }
    cached
  }
})

# what bayesplot itself does with a km type, given only y and yrep,
# which is all brms and frmtmb pass it
ppt_km_verdict <- function() {
  tryCatch({
    p <- bayesplot::ppc_km_overlay(c(1, 2, 3),
                                   matrix(c(1, 2, 3, 2, 3, 1), nrow = 2,
                                          byrow = TRUE))
    invisible(ggplot2::ggplot_build(p))
    "ok"
  }, error = function(e) "bayesplot")
}

# What brms does with each type on each response, as measured. "ok" is
# a plot; "bayesplot" is bayesplot's own refusal of the data, which
# brms passes through and so does this package; "loo" is a type that
# weights posterior draws, which brms answers and a maximum-likelihood
# fit refuses with the reason. km_* need ggfortify, which brms and
# bayesplot both require.
ppt_verdict <- function(type, resp) {
  if (startsWith(type, "loo_")) return("loo")
  if (startsWith(type, "km_")) {
    # the km types need a package bayesplot suggests, and a status
    # argument neither brms nor frmtmb passes. Which of the two
    # bayesplot complains about first is this machine's business, so
    # the verdict is READ OFF bayesplot with the arguments both
    # packages pass, rather than from a skip
    return(ppt_km_verdict())
  }
  discrete <- resp %in% c("cnt", "bin")
  if (type %in% c("bars", "bars_grouped", "rootogram",
                  "rootogram_grouped")) {
    return(if (discrete) "ok" else "bayesplot")
  }
  if (type %in% c("calibration", "calibration_grouped")) {
    return(if (identical(resp, "bin")) "ok" else "bayesplot")
  }
  if (type %in% c("calibration_overlay", "calibration_overlay_grouped")) {
    # brms passes no `prep`, which these two require
    return("bayesplot")
  }
  "ok"
}

test_that("every bayesplot ppc type does on a fit what it does in brms", {
  d <- ppt_data()
  fits <- list(
    y = frm(bf(y ~ x + (1 | g)) + gaussian(), data = d),
    cnt = frm(bf(cnt ~ x + (1 | g)) + poisson(), data = d),
    bin = frm(bf(bin ~ x + (1 | g)) + bernoulli(), data = d)
  )
  types <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  # every type bayesplot lists, not a sample: 49 in bayesplot 1.16.0
  expect_gte(length(types), 49L)
  ns <- asNamespace("bayesplot")
  seen <- character(0)
  for (resp in names(fits)) {
    for (ty in types) {
      fa <- names(formals(get(paste0("ppc_", ty), envir = ns)))
      args <- list(fits[[resp]], type = ty, ndraws = 10)
      if ("group" %in% fa) args$group <- "g"
      if ("x" %in% fa) args$x <- "x"
      want <- ppt_verdict(ty, resp)
      set.seed(1)
      got <- tryCatch({
        p <- suppressWarnings(suppressMessages(do.call(pp_check, args)))
        invisible(suppressWarnings(suppressMessages(
          ggplot2::ggplot_build(p))))
        "ok"
      }, frmtmb_error = function(e) {
        if (grepl("Pareto-smoothed", conditionMessage(e))) "loo" else
          paste("frmtmb:", conditionMessage(e))
      }, error = function(e) "bayesplot")
      expect_identical(got, want, info = paste(resp, ty))
      seen <- c(seen, paste(resp, ty))
    }
  }
  # the loop really ran every pair
  expect_identical(length(seen), length(types) * 3L)
})

# Three blocks, not one, for the same reason as the refusals above: on
# 0.61.0 the first pp_check() call here ERRORS, so in one block the
# other two types were never reached on the base arm.
test_that("a grouped type panels by the model's own variable", {
  d <- ppt_data()
  fit <- ppt_gauss_fit()
  set.seed(2)
  p <- pp_check(fit, type = "stat_grouped", group = "g", ndraws = 10)
  # one panel per level, and the rows of each panel are that level's
  expect_setequal(as.character(unique(p$data$group)), levels(d$g))
})

test_that("a grouped type keeps the observations in row order", {
  d <- ppt_data()
  fit <- ppt_gauss_fit()
  set.seed(2)
  pv <- pp_check(fit, type = "violin_grouped", group = "g", ndraws = 10)
  yobs <- pv$data[pv$data$is_y, ]
  yobs <- yobs[order(yobs$y_id), ]
  expect_identical(nrow(yobs), nrow(d))
  expect_identical(as.character(yobs$group), as.character(d$g))
})

test_that("an x variable arrives as the model's values, in row order", {
  d <- ppt_data()
  fit <- ppt_gauss_fit()
  set.seed(2)
  pi <- pp_check(fit, type = "intervals", x = "x", ndraws = 10)
  expect_equal(pi$data$x, d$x)
})

# One refusal per block. They were one block until the recheck of punch
# round 1, and on 0.61.0 the first of them ERRORS rather than failing
# ("object 'ppc_violin' not found" is not a frmtmb_error), so testthat
# ended the block there and the other eight never ran on the base arm.
# A count from such a block says nothing about the assertions behind the
# error, which is why each one now stands alone.
ppt_refusals <- list(
  list(what = "a type bayesplot does not have",
       # bayesplot has ppc_violin_grouped and no ppc_violin
       call = quote(pp_check(fit, type = "violin")),
       pats = c("'violin' is not a valid ppc type", "'violin_grouped'")),
  list(what = "a grouped type called without a group",
       call = quote(pp_check(fit, type = "stat_grouped")),
       pats = "Argument 'group' is required"),
  list(what = "a group name the data does not carry",
       call = quote(pp_check(fit, type = "stat_grouped",
                             group = "nosuch")),
       pats = "'nosuch' could not be found"),
  list(what = "a group naming a variable the model does not use",
       # z is in the data and not in the model, which brms refuses too
       call = quote(pp_check(fit, type = "stat_grouped", group = "z")),
       pats = "'z' could not be found"),
  list(what = "an x naming a variable the model does not use",
       call = quote(pp_check(fit, type = "intervals", x = "z")),
       pats = "'z' could not be found"),
  list(what = "a group that is not one string",
       call = quote(pp_check(fit, type = "stat_grouped",
                             group = c("g", "g"))),
       pats = "single string"),
  list(what = "a newdata that lacks the response, for a ppc type",
       # brms refuses it too: "Response variables must be specified in
       # 'newdata'" (dev/simnewdata-log/brms.txt)
       call = quote(pp_check(fit, type = "dens_overlay",
                             newdata = ppt_data()[1:10, c("x", "g")])),
       pats = "does not supply 'y'"),
  list(what = "draw_ids, which a maximum-likelihood fit has none of",
       call = quote(pp_check(fit, ndraws = 5, draw_ids = 1:3)),
       pats = "draw_ids")
)

for (ppt_i in seq_along(ppt_refusals)) {
  local({
    case <- ppt_refusals[[ppt_i]]
    test_that(paste("pp_check() refuses", case$what), {
      fit <- ppt_gauss_fit()
      for (pat in case$pats) {
        expect_error(eval(case$call), pat, class = "frmtmb_error",
                     info = pat)
      }
    })
  })
}

# These four live in a block of their OWN on purpose. In the block above
# they never ran against 0.61.0: its `type = "violin"` assertion ERRORS
# there, testthat ends the block, and every assertion after it is
# skipped, so a count from that block says nothing about these. Each one
# fails on 0.61.0 when it is reached (dev/correct-log/punch1-base.txt).
test_that("pp_check() passes `resp` and an unused `group` as brms does", {
  d <- ppt_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
  # `resp` on a one-response model is accepted and ignored, whatever it
  # names: brms's validate_resp() returns NULL there, so brms plots
  set.seed(3)
  p1 <- suppressMessages(pp_check(fit, type = "dens_overlay", ndraws = 5))
  set.seed(3)
  # SILENTLY: `resp` is a formal here, so it does not fall into the dots
  # and make bayesplot warn, which is what brms does and what 0.61.0
  # did not
  p2 <- expect_no_warning(suppressMessages(
    pp_check(fit, type = "dens_overlay", ndraws = 5, resp = "w")))
  expect_equal(p1$data, p2$data)
  # a type with no `group` argument gets the COLUMN in its dots and
  # bayesplot warns, which is what brms does with the same call
  expect_warning(pp_check(fit, type = "dens_overlay", group = "g",
                          ndraws = 5),
                 "unrecognized and ignored")
  # ... and a group name the model does not use is not looked up for
  # such a type, so nothing is passed and nothing warns, again as in brms
  expect_silent(suppressMessages(
    pp_check(fit, type = "dens_overlay", group = "nosuch", ndraws = 5)))
  # `x` is the ASYMMETRIC half. brms writes `as.numeric(data[[x]])`
  # where it writes a bare `data[[group]]`, so an x name the data does
  # not carry becomes numeric(0) rather than NULL, stays in the argument
  # list and WARNS, where the same group name is silent. Measured on
  # brms in dev/correct-log/punch3-brms-xarg.txt; both halves are
  # asserted because copying one of them was the defect.
  expect_warning(pp_check(fit, type = "dens_overlay", x = "x",
                          ndraws = 5),
                 "unrecognized and ignored")
  expect_warning(pp_check(fit, type = "dens_overlay", x = "nosuch",
                          ndraws = 5),
                 "unrecognized and ignored")
})

test_that("prefix = \"ppd\" plots the simulated responses alone", {
  d <- ppt_data()
  fit <- frm(bf(cnt ~ x + (1 | g)) + poisson(), data = d)
  for (ty in c("dens_overlay", "stat_grouped", "intervals")) {
    fa <- names(formals(get(paste0("ppd_", ty), asNamespace("bayesplot"))))
    args <- list(fit, type = ty, ndraws = 10, prefix = "ppd")
    if ("group" %in% fa) args$group <- "g"
    if ("x" %in% fa) args$x <- "x"
    p <- suppressMessages(do.call(pp_check, args))
    expect_s3_class(p, "ggplot")
    # no observed series on a ppd plot
    pd <- as.data.frame(p$data)
    expect_false(any(c("is_y", "y_obs") %in% names(pd)), info = ty)
    expect_false(any(pd[["variable"]] %in% "y"), info = ty)
  }
  expect_error(pp_check(fit, type = "violin_grouped", prefix = "ppd",
                        group = "g"),
               "not a valid ppd type", class = "frmtmb_error")
})

test_that("error_binned is refused on a category response, as in brms", {
  set.seed(9)
  n <- 150
  d <- data.frame(x = stats::rnorm(n))
  d$o <- factor(cut(0.8 * d$x + stats::rlogis(n), c(-Inf, -0.5, 0.7, Inf),
                    labels = FALSE), ordered = TRUE)
  d$k <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  fo <- frm(bf(o ~ x) + cumulative(), data = d)
  fk <- frm(bf(k ~ x) + categorical(), data = d)
  for (f in list(fo, fk)) {
    expect_error(pp_check(f, type = "error_binned", ndraws = 5),
                 "not available for polytomous", class = "frmtmb_error")
    # brms answers error_hist on both, on the category codes
    expect_s3_class(pp_check(f, type = "error_hist", ndraws = 5), "ggplot")
  }
  # the categorical draws reach bayesplot as the 1..K codes y is on
  set.seed(3)
  p <- pp_check(fk, type = "bars", ndraws = 5)
  expect_setequal(p$data$x, 1:3)
})

# pp_check(newdata = ) answers as brms answers, rather than refusing
# (lane wt-simnewdata). brms 2.23.0 on its own example fit: a ggplot
# whose observed series is the newdata's 10 responses, and a grouped
# type whose groups are the newdata's (dev/simnewdata-log/brms.txt).
test_that("pp_check(newdata = ) plots the newdata's rows", {
  fit <- ppt_gauss_fit()
  nd <- ppt_data()[1:10, ]
  set.seed(4)
  p <- pp_check(fit, newdata = nd, ndraws = 5)
  expect_s3_class(p, "ggplot")
  pd <- as.data.frame(p$data)
  expect_equal(pd$value[pd$is_y], nd$y)
  expect_identical(nrow(pd), 10L * 6L)
})

test_that("pp_check(newdata = ) reads `group` off the newdata", {
  fit <- ppt_gauss_fit()
  nd <- ppt_data()[1:12, ]
  nd$g <- factor(nd$g, levels = levels(nd$g))
  nd <- nd[nd$g %in% c("a", "b"), ]
  set.seed(5)
  p <- pp_check(fit, type = "violin_grouped", group = "g", newdata = nd,
                ndraws = 5)
  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(unique(p$data$group)), c("a", "b"))
})

test_that("pp_check(newdata = , prefix = \"ppd\") needs no response", {
  fit <- ppt_gauss_fit()
  nd <- ppt_data()[1:10, c("x", "g")]
  set.seed(6)
  p <- pp_check(fit, newdata = nd, ndraws = 5, prefix = "ppd")
  expect_s3_class(p, "ggplot")
})

test_that("pp_check(newdata = ) drops a missing response as brms does", {
  fit <- ppt_gauss_fit()
  nd <- ppt_data()[1:10, ]
  nd$y[3] <- NA
  set.seed(7)
  expect_warning(p <- pp_check(fit, newdata = nd, ndraws = 5),
                 "NA responses are not shown")
  pd <- as.data.frame(p$data)
  expect_equal(pd$value[pd$is_y], nd$y[-3])
})

test_that("pp_check(re_formula = ~1) redraws the group effects, as NA does", {
  # ~1 used to condition on the fitted effects in simulate() while
  # meaning "no group-level effects" in predict()
  fit <- ppt_gauss_fit()
  set.seed(8)
  p1 <- pp_check(fit, ndraws = 5, re_formula = ~1)
  set.seed(8)
  pna <- pp_check(fit, ndraws = 5, re_formula = NA)
  expect_identical(p1$data, pna$data)
})

test_that("error_binned on a multinomial fit calls it counts, not a category", {
  set.seed(105)
  n <- 80
  dm <- data.frame(x = stats::rnorm(n), n = 15L)
  dm$Y <- t(vapply(seq_len(n), function(i) {
    as.vector(stats::rmultinom(1, 15, c(0.3, 0.3, 0.4)))
  }, numeric(3)))
  fit <- frm(bf(Y | trials(n) ~ x) + multinomial(K = 3), data = dm)
  expect_error(pp_check(fit, type = "error_binned"),
               "response is a set of counts over categories",
               class = "frmtmb_error")
})
