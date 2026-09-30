# cmc, brms's switch for cell-mean coding: with cmc = FALSE a formula
# without an intercept keeps a factor's treatment contrasts and drops
# only the intercept column (brms:::validate_terms()). Checked against
# brms 2.23.0's standata(), which needs no Stan compile.

cmc_strip <- function(X) {
  X <- as.matrix(X)
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  unname(X)
}

test_that("bf(cmc = FALSE) builds brms's population-level design", {
  skip_on_cran()
  skip_if_not_installed("brms")
  df <- data.frame(y = 1:10, g = rep(c("a", "b"), 5))
  # brms's own case, whose (1 | y) is one level per row
  fr <- allow_warnings(
    frm(bf(y ~ 0 + g + (1 | y), cmc = FALSE), data = df, dry_run = "frame"),
    "gives every observation its own random effect")
  sd <- brms::standata(brms::bf(y ~ 0 + g + (1 | y), cmc = FALSE), df)
  X <- fr$linpreds[["y.mu"]]$X
  expect_identical(colnames(X), colnames(sd$X))
  expect_equal(cmc_strip(X), unname(sd$X), ignore_attr = TRUE)
  # R's rule without cmc = FALSE: one column per level
  fr0 <- frm(bf(y ~ 0 + g), data = df, dry_run = "frame")
  expect_identical(colnames(fr0$linpreds[["y.mu"]]$X), c("ga", "gb"))
  # a formula with an intercept does not change
  fr1 <- frm(bf(y ~ g, cmc = FALSE), data = df, dry_run = "frame")
  expect_identical(colnames(fr1$linpreds[["y.mu"]]$X),
                   c("(Intercept)", "gb"))
})

test_that("lf(cmc = FALSE) reaches a parameter's group-level terms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(5)
  d3 <- data.frame(g = factor(rep(c("a", "b", "c"), 20)),
                   h = factor(rep(1:6, each = 10)))
  d3$y <- rnorm(60, as.numeric(d3$g))
  fr <- allow_warnings(
    frm(bf(y ~ 1) + lf(sigma ~ 0 + g + (0 + g | h), cmc = FALSE),
        data = d3, dry_run = "frame"),
    "confounded")
  sd <- brms::standata(brms::bf(y ~ 1) +
                         brms::lf(sigma ~ 0 + g + (0 + g | h), cmc = FALSE),
                       d3)
  lp <- fr$linpreds[["y.sigma"]]
  expect_identical(colnames(lp$X), colnames(sd$X_sigma))
  expect_equal(cmc_strip(lp$X), unname(sd$X_sigma), ignore_attr = TRUE)
  cp <- fr$re_blocks[[1]]$components[[1]]
  expect_identical(cp$cnms, c("gb", "gc"))
  # brms has two coefficients per level, Z_1_sigma_1 and Z_1_sigma_2,
  # and no third; the frame's Z is level-major
  expect_false("Z_1_sigma_3" %in% names(sd))
  Z <- as.matrix(lp$Z)
  expect_equal(rowSums(Z[, seq(1, ncol(Z), by = 2)]),
               as.numeric(sd$Z_1_sigma_1), ignore_attr = TRUE)
  expect_equal(rowSums(Z[, seq(2, ncol(Z), by = 2)]),
               as.numeric(sd$Z_1_sigma_2), ignore_attr = TRUE)
})

test_that("a cmc = FALSE fit is the model of its design", {
  set.seed(5)
  d3 <- data.frame(g = factor(rep(c("a", "b", "c"), 80)),
                   h = factor(rep(1:16, each = 15)))
  u <- matrix(rnorm(32), 16, 2)
  d3$y <- rnorm(240, as.numeric(d3$g), 0.5) +
    (d3$g == "b") * u[d3$h, 1] + (d3$g == "c") * u[d3$h, 2]
  fit <- frm(bf(y ~ 0 + g, cmc = FALSE), data = d3)
  # the same columns written out by hand
  d3$gb <- as.numeric(d3$g == "b")
  d3$gc <- as.numeric(d3$g == "c")
  hand <- frm(y ~ 0 + gb + gc, data = d3)
  expect_identical(names(fixef_by_dpar(fit)$mu), c("gb", "gc"))
  expect_identical(as.numeric(logLik(fit)), as.numeric(logLik(hand)))
  # and a prediction on new data rebuilds the same columns, group-level
  # terms included
  fit_re <- frm(bf(y ~ 0 + g + (0 + g | h), cmc = FALSE), data = d3)
  expect_equal(frm_linpred(fit_re, newdata = d3), frm_linpred(fit_re),
               ignore_attr = TRUE)
})

test_that("cmc = FALSE reaches a multi-membership term, as in brms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(21)
  n <- 120
  d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                  h = factor(sample(1:8, n, TRUE)),
                  h2 = factor(sample(1:8, n, TRUE)))
  d$y <- rnorm(n)
  fr <- frm(bf(y ~ x + (0 + f | mm(h, h2)), cmc = FALSE), data = d,
            dry_run = "frame")
  cp <- fr$re_blocks[[1]]$components[[1]]
  expect_identical(cp$cnms, c("fb", "fc"))
  sd <- brms::standata(brms::bf(y ~ x + (0 + f | mm(h, h2)), cmc = FALSE),
                       d)
  # brms: two coefficients per member, each weighted by W_1_<member>
  expect_false("Z_1_3_1" %in% names(sd))
  Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  for (k in 1:2) {
    b <- sd$W_1_1 * sd[[paste0("Z_1_", k, "_1")]] +
      sd$W_1_2 * sd[[paste0("Z_1_", k, "_2")]]
    expect_equal(rowSums(Z[, seq(k, ncol(Z), by = 2)]), as.numeric(b),
                 ignore_attr = TRUE)
  }
  # a prediction rebuilds the members' designs the same way
  fit <- allow_warnings(
    frm(bf(y ~ x + (0 + f | mm(h, h2)), cmc = FALSE), data = d),
    c("convergence", "gradient", "singular"))
  expect_equal(frm_linpred(fit, newdata = d), frm_linpred(fit),
               ignore_attr = TRUE)
})

test_that("cmc = FALSE is refused on a structure read level by level", {
  set.seed(22)
  d <- data.frame(t = factor(rep(1:4, 15)), g = factor(rep(1:15, each = 4)))
  d$y <- rnorm(60)
  expect_error(frm(bf(y ~ 1 + ar1(0 + t | g), cmc = FALSE), data = d,
                   dry_run = "frame"),
               "cmc = FALSE does not apply to ar1", class = "frmtmb_error")
})

test_that("cmc = FALSE is accepted where it changes no column", {
  # numeric slopes, and an interaction without its factor's main effect,
  # have the same columns under either setting, so the structures that
  # read one coefficient per level have nothing to lose
  # (dev/reviews/2026-09-29-formula2.md, B3)
  set.seed(41)
  d <- data.frame(x1 = rnorm(60), x2 = rnorm(60),
                  f = factor(rep(c("a", "b", "c"), 20)),
                  g = factor(rep(1:15, each = 4)))
  d$y <- rnorm(60)
  cnms <- function(f) {
    fr <- frm(f, data = d, dry_run = "frame")
    fr$re_blocks[[1]]$cnms
  }
  for (cs in c("cs", "ar1", "hetar1", "homcs", "toep", "homtoep")) {
    t_on <- str2lang(paste0(cs, "(0 + x1 + x2 | g)"))
    f_true <- stats::as.formula(call("~", quote(y), call("+", 1, t_on)))
    f_false <- bf(f_true, cmc = FALSE)
    expect_identical(cnms(f_false), cnms(f_true), label = cs)
    expect_identical(cnms(f_false), c("x1", "x2"), label = cs)
  }
  expect_identical(cnms(bf(y ~ 1 + cs(0 + f:x1 | g), cmc = FALSE)),
                   cnms(y ~ 1 + cs(0 + f:x1 | g)))
})

test_that("re_formula names a cmc = FALSE term as it was fitted", {
  set.seed(5)
  d3 <- data.frame(f = factor(rep(c("a", "b", "c"), 80)),
                   h = factor(rep(1:16, each = 15)))
  u <- matrix(rnorm(32), 16, 2)
  d3$y <- rnorm(240, as.numeric(d3$f), 0.5) +
    (d3$f == "b") * u[d3$h, 1] + (d3$f == "c") * u[d3$h, 2]
  fit <- frm(bf(y ~ 0 + f + (0 + f | h), cmc = FALSE), data = d3)
  # the only group-level term, so naming it is the full prediction
  expect_equal(fitted(fit, re_formula = ~ (0 + f | h)), fitted(fit),
               ignore_attr = TRUE)
})

test_that("print() shows cmc", {
  expect_output(print(bf(y ~ 0 + g, cmc = FALSE)),
                "y ~ 0 + g (cmc = FALSE)", fixed = TRUE)
  expect_output(print(lf(sigma ~ 0 + g, cmc = FALSE)),
                "sigma ~ 0 + g (cmc = FALSE)", fixed = TRUE)
})

test_that("cmc is a flag, and lf() says what it takes", {
  expect_error(bf(y ~ 0 + g, cmc = "no"), "`cmc` must be TRUE or FALSE",
               class = "frmtmb_error")
  expect_error(lf(sigma ~ 0 + g, cmc = NA), "`cmc` must be TRUE or FALSE",
               class = "frmtmb_error")
  expect_identical(attr(lf(sigma ~ 0 + g, cmc = FALSE)$pforms$sigma, "cmc"),
                   FALSE)
  # the message a non-formula gets is true of lf() now
  expect_error(lf(sigma ~ x, 3), "lf[(][)] takes parameter formulas",
               class = "frmtmb_error")
  # bf() on a formula already built takes cmc as well
  expect_false(bf(bf(y ~ 0 + g), cmc = FALSE)$cmc)
  # a misspelled switch is still refused, not taken as a dpar setting
  # (dev/argspell-findings.md): a flag is no constant, and a string is
  # an equation only between two names of one class
  expect_error(bf(y ~ 0 + g, cmcc = FALSE),
               "Cannot interpret bf[(][)] argument 'cmcc'",
               class = "frmtmb_error")
  expect_error(lf(sigma ~ 0 + g, cmcc = FALSE),
               "lf[(][)] takes parameter formulas", class = "frmtmb_error")
  expect_error(bf(y ~ x, famly = "gaussian"),
               "has no argument `famly`. Did you mean `family`?",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(lf(sigma ~ x, rsp = "y1"),
               "has no argument `rsp`. Did you mean `resp`?",
               fixed = TRUE, class = "frmtmb_error")
})
