# conditional_effects() against brms 2.23.0 where a call names effects,
# asks an ordinal fit for the expected category number, or displays a
# model with mi(x, idx = ). brms's answers are dev/ceplot-log/brms-mi.txt
# (dev/ceplot-brms-mi.R); 0.66.0's are dev/ceplot-log/probe-base.txt.

parity_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(1)
      dd <- data.frame(x = stats::rnorm(100),
                       f = factor(rep(c("a", "b"), 50)),
                       z = stats::rnorm(100))
      dd$y <- stats::rnorm(100, 1 + 0.5 * dd$x + (dd$f == "b"))
      cache <<- frm(bf(y ~ x * f), family = gaussian(), data = dd)
    }
    cache
  }
})

test_that("one invalid effect among valid ones is dropped with a warning", {
  # brmsfit-methods:205. 0.66.0 stopped: "Variable 'xx' is not stored in
  # the model frame"
  fit <- parity_fit()
  ce <- allow_warnings(
    conditional_effects(fit, effects = c("xx", "x")),
    "Some specified effects are invalid for this model",
    require = "Some specified effects are invalid for this model: 'xx'")
  expect_named(ce, "x")
  # a variable in the data that the model does not read is invalid too,
  # as in brms: its curve would be flat by construction
  ce <- allow_warnings(
    conditional_effects(fit, effects = c("z", "f")),
    "Some specified effects are invalid for this model",
    require = "invalid for this model: 'z'")
  expect_named(ce, "f")
  # a pair with an invalid side, and a pair naming one variable twice
  ce <- allow_warnings(
    conditional_effects(fit, effects = c("x:xx", "x:x", "f")),
    "Some specified effects are invalid for this model",
    require = "invalid for this model: 'x:xx', 'x:x'")
  expect_named(ce, "f")
  # all invalid: brms's error, with the valid variables listed
  expect_error(conditional_effects(fit, effects = "xx"),
               "All specified effects are invalid for this model.")
  expect_error(conditional_effects(fit, effects = c("xx", "z")),
               "Valid effects are (combinations of): 'x', 'f'",
               fixed = TRUE)
  # the order-three check comes first, as in brms (brmsfit-methods:207)
  expect_error(conditional_effects(fit, effects = "x:f:xx"),
               "at most two variables")
  # the guard absent: valid effects, including a grouping variable of a
  # group-level term and a pair, pass with no warning
  set.seed(2)
  d <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- stats::rnorm(60, d$x + stats::rnorm(6)[d$g])
  fg <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  expect_no_warning(ce <- conditional_effects(fg, effects = c("x", "g",
                                                               "x:g"),
                                              resolution = 3))
  expect_named(ce, c("x", "g", "x:g"))
})

test_that("the expected category number of an ordinal fit warns as brms", {
  # brmsfit-methods:217. brms warns on its default call, whose display is
  # this one; frmtmb's default is the per-category display and warns
  # there only when categorical = FALSE asks for the summary
  set.seed(2)
  do <- data.frame(x = stats::rnorm(200))
  do$y <- factor(cut(do$x + stats::rlogis(200), c(-Inf, -1, 0, 1, Inf)),
                 ordered = TRUE)
  fo <- frm(bf(y ~ x), family = cumulative(), data = do)
  msg <- "Predictions are treated as continuous variables"
  ce <- allow_warnings(conditional_effects(fo, categorical = FALSE,
                                           resolution = 5),
                       msg, require = msg)
  expect_named(ce, "x")
  # the guard absent: the per-category display and the latent predictor
  # treat nothing as continuous, and brms gives neither a warning
  expect_no_warning(conditional_effects(fo, resolution = 5))
  expect_no_warning(conditional_effects(fo, dpar = "mu", resolution = 5))
})

test_that("a display of mi(x, idx = ) is refused by name, with the reason", {
  # 0.66.0 stopped inside the prediction with "Could not match all
  # indices"; brms stops with "Index of response 'x' contains duplicated
  # values." (dev/ceplot-log/brms-mi.txt)
  set.seed(26)
  n <- 120
  dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                   s = rep(c(TRUE, FALSE), n / 2), w = stats::rnorm(n))
  dm$x <- stats::rnorm(n)
  dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w +
    stats::rnorm(n, sd = 0.5)
  mv <- frm(bf(y ~ mi(x, idx = g1) + w) +
              bf(x | mi() + index(g2) + subset(s) ~ 1),
            data = dm, family = gaussian())
  expect_error(conditional_effects(mv, "w", resp = "y"),
               "cannot display response 'y': its predictor reads mi\\(x, idx")
  expect_error(conditional_effects(mv, resp = "y"), "the same row")
  # the guard absent: a response that does not read the term answers,
  # and so does mi() without idx
  expect_named(conditional_effects(mv, "w", resp = "x", resolution = 3),
               "w")
  set.seed(3)
  dn <- data.frame(w = stats::rnorm(80))
  dn$x <- stats::rnorm(80, 0.5 * dn$w)
  dn$y <- stats::rnorm(80, 0.5 * dn$x + 0.3 * dn$w)
  dn$x[c(3, 9)] <- NA
  mn <- frm(bf(y ~ mi(x) + w) + bf(x | mi() ~ w), data = dn,
            family = gaussian())
  expect_named(conditional_effects(mn, "w", resp = "y", resolution = 3),
               "w")
})

# Punch round 1.

test_that("a cs() variable is a valid effect", {
  # review m2: the valid set without the cs() variables survived this
  # file; brms's get_all_effects() counts a cs() term's variables
  set.seed(405)
  n <- 300
  d <- data.frame(x = stats::rnorm(n),
                  fc = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$yo <- factor(cut(d$x + c(a = -1, b = 0, c = 1.5)[as.character(d$fc)] +
                       stats::rlogis(n), c(-Inf, -0.5, 0.8, Inf),
                     labels = FALSE), ordered = TRUE)
  ff <- frm(bf(yo ~ x + cs(fc)), family = sratio(), data = d)
  expect_no_warning(ce <- conditional_effects(ff, effects = c("fc", "x"),
                                              resolution = 3))
  expect_named(ce, c("fc:cats__", "x:cats__"))
})

test_that("a nonlinear parameter's name is not a valid effect", {
  # review m3: `a` and `b` were listed as valid, and effects = c("a",
  # "x") stopped with "Variable 'a' is not stored"; brms drops `a` with
  # its warning and lists x and z
  set.seed(11)
  d <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
  d$yp <- 2 * exp(0.3 * d$x) + 0.2 * d$z + stats::rnorm(150, 0, 0.3)
  fn <- frm(bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
            family = gaussian(), data = d)
  seen <- character()
  ce <- withCallingHandlers(
    conditional_effects(fn, effects = c("a", "x"), resolution = 3,
                        band = "boot", boot = 5, seed = 1),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_named(ce, "x")
  expect_length(seen, 1L)
  expect_match(seen, "invalid for this model: 'a'", fixed = TRUE)
  valid <- sub(".*combinations of[)]: ", "", seen)
  expect_setequal(strsplit(valid, ", ", fixed = TRUE)[[1L]],
                  c("'x'", "'z'"))
})

test_that("method = \"predict\" on an ordinal fit stops before it warns", {
  # review m8: the ordinal warning came first and the call then stopped
  set.seed(2)
  do <- data.frame(x = stats::rnorm(200))
  do$y <- factor(cut(do$x + stats::rlogis(200), c(-Inf, -1, 0, 1, Inf)),
                 ordered = TRUE)
  fo <- frm(bf(y ~ x), family = cumulative(), data = do)
  expect_no_warning(expect_error(
    conditional_effects(fo, categorical = FALSE, method = "predict"),
    "has no meaning on an ordinal family"))
})

test_that("a variable read through poly() or a transform is an effect", {
  # pre-existing: the model frame keeps poly(x, 2) and log(abs(z) + 1)
  # but not x and z, so naming either stopped with "not stored in the
  # model frame" and the default display found nothing; brms accepts
  # both
  set.seed(11)
  d <- data.frame(x = stats::rnorm(80), z = stats::rnorm(80))
  d$y <- stats::rnorm(80, 1 + d$x - 0.5 * d$x^2 + log(abs(d$z) + 1))
  fit <- frm(bf(y ~ poly(x, 2) + log(abs(z) + 1)), family = gaussian(),
             data = d)
  expect_named(conditional_effects(fit, resolution = 5), c("x", "z"))
  ce <- conditional_effects(fit, c("x", "z"), resolution = 5)
  # the curve is the fit's at the grid: x varies with z at its mean
  nd <- data.frame(x = ce$x$x, z = mean(d$z))
  f <- fitted(fit, newdata = nd)[, "Estimate"]
  eps <- .Machine$double.eps
  expect_lt(max(abs(ce$x$estimate__ - f)), 64 * eps * max(abs(f)))
  expect_equal(range(ce$z$z), range(d$z))
})
