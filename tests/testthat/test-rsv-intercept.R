# brms's reserved `Intercept` (`y ~ 0 + Intercept + x`) and
# bf(center = FALSE). The likelihood is that of `y ~ 1 + x`; what brms
# changes is the prior class of the intercept ("b", not "Intercept"),
# the centering a class "Intercept" prior implies, and the coefficient
# order. dev/icpt0-findings.md has the validation, and
# dev/icpt0-validate.R the script behind its numbers.

ri_data <- local({
  set.seed(1)
  n <- 80
  d <- data.frame(x = stats::rnorm(n, 3), z = stats::rnorm(n),
                  f = factor(rep(letters[1:4], n / 4)))
  d$y <- 1 + 0.5 * d$x + 0.4 * as.numeric(d$f) +
    stats::rnorm(n, sd = exp(0.2 * d$z))
  d$cnt <- stats::rpois(n, exp(0.2 + 0.2 * d$x))
  d
})

ri_x <- function(form, dpar = "mu") {
  fr <- frm(form, data = ri_data, dry_run = "frame")
  for (lp in fr$linpreds) if (identical(lp$dpar, dpar)) return(lp$X)
}

ri_rows <- function(tab) {
  tab <- as.data.frame(tab)
  tab <- tab[tab$class %in% c("b", "Intercept"), ]
  sort(paste(tab$class, tab$coef, tab$dpar, tab$nlpar, tab$resp, sep = "|"))
}

test_that("0 + Intercept is the model 1 + x under ML, bitwise", {
  # an identity: the rewritten formula builds the same design, the same
  # start and the same tape, so the optimizer takes the same path
  f0 <- frm(bf(y ~ 0 + Intercept + x + f), data = ri_data)
  f1 <- frm(bf(y ~ x + f), data = ri_data)
  expect_identical(logLik(f0), logLik(f1))
  expect_identical(fixef(f0), fixef(f1))
  s0 <- frm(bf(y ~ x, sigma ~ 0 + Intercept + z), data = ri_data)
  s1 <- frm(bf(y ~ x, sigma ~ z), data = ri_data)
  # the same numbers, where brms's order puts sigma_Intercept (see below)
  expect_identical(fixef(s0), fixef(s1)[rownames(fixef(s0)), ])
  # written last, the column moves and the path with it, so equality is
  # to the optimizer's precision, judged against the standard errors
  f2 <- frm(bf(y ~ 0 + x + f + Intercept), data = ri_data)
  fe <- fixef(f1)
  expect_lt(abs(as.numeric(logLik(f2) - logLik(f1))),
            1e-10 * abs(as.numeric(logLik(f1))))
  expect_lt(max(abs(fixef(f2)[rownames(fe), 1] - fe[, 1]) / fe[, 2]), 1e-4)
})

test_that("a data column of ones named Intercept gives brms's model", {
  # The pin. Before, `Intercept` was an ordinary covariate, so this data
  # column fitted, but R's no-intercept rule coded f with all four
  # levels beside it: the design was aliased, a column was dropped, and
  # the coefficients were not brms's. Prediction needed the column too.
  d <- ri_data
  d$Intercept <- 1
  fit <- frm(bf(y ~ 0 + Intercept + f), data = d)
  ref <- frm(bf(y ~ f), data = d)
  expect_identical(fixef(fit), fixef(ref))
  nd <- data.frame(f = factor(c("b", "d"), levels = letters[1:4]))
  expect_identical(fitted(fit, newdata = nd), fitted(ref, newdata = nd))
})

test_that("the design and the prior rows are brms's", {
  skip_on_cran()
  skip_if_not_installed("brms")
  forms <- list(
    y ~ 0 + Intercept + x, y ~ 0 + x + Intercept, y ~ x + Intercept - 1,
    y ~ 0 + Intercept + f, y ~ 0 + x * f + Intercept,
    y ~ 0 + Intercept + x:f)
  for (fm in forms) {
    xb <- brms::make_standata(fm, data = ri_data)$X
    xf <- ri_x(bf(fm))
    expect_identical(sub("^[(]Intercept[)]$", "Intercept", colnames(xf)),
                     colnames(xb), label = deparse1(fm))
    expect_identical(as.vector(as.matrix(xf)), as.vector(unclass(xb)),
                     label = deparse1(fm))
  }
  expect_identical(
    ri_rows(get_prior(bf(y ~ x, sigma ~ 0 + Intercept + z),
                      data = ri_data)),
    ri_rows(brms::get_prior(brms::bf(y ~ x, sigma ~ 0 + Intercept + z),
                            data = ri_data)))
  expect_identical(
    ri_rows(get_prior(bf(y ~ 0 + Intercept + f), data = ri_data)),
    ri_rows(brms::get_prior(y ~ 0 + Intercept + f, data = ri_data)))
  expect_identical(
    ri_rows(get_prior(bf(y ~ x, center = FALSE) +
                        lf(sigma ~ z, center = FALSE), data = ri_data)),
    ri_rows(brms::get_prior(brms::bf(y ~ x, center = FALSE) +
                              brms::lf(sigma ~ z, center = FALSE),
                            data = ri_data)))
})

test_that("factors take treatment contrasts beside the reserved Intercept", {
  # R's own rule gives the first factor all its levels in a formula
  # without an intercept, which next to a column of ones is aliased;
  # brms fits the intercept design and so does frmtmb
  x <- ri_x(bf(y ~ 0 + Intercept + f))
  expect_identical(colnames(x), c("(Intercept)", "fb", "fc", "fd"))
  expect_identical(colnames(ri_x(bf(y ~ 0 + x * f + Intercept))),
                   c("x", "fb", "fc", "fd", "(Intercept)", "x:fb", "x:fc",
                     "x:fd"))
})

test_that("a class b prior reaches the intercept; class Intercept does not", {
  s <- 0.1
  fit <- frm(bf(y ~ 0 + Intercept + x), data = ri_data,
             prior = set_prior("normal(0, 0.1)", class = "b"))
  # brms: lprior += normal_lpdf(b | 0, 0.1) over (b_Intercept, b_x), and
  # the likelihood of y ~ 1 + x; its mode, written out independently
  X <- cbind(1, ri_data$x)
  nlp <- function(th) {
    -(sum(stats::dnorm(ri_data$y, drop(X %*% th[1:2]), exp(th[3]),
                       log = TRUE)) +
        sum(stats::dnorm(th[1:2], 0, s, log = TRUE)))
  }
  o <- stats::optim(c(0, 0, 0), nlp, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 1000))
  fe <- fixef(fit)
  expect_lt(max(abs(fe[c("Intercept", "x"), 1] - o$par[1:2]) /
                  fe[c("Intercept", "x"), 2]), 1e-3)
  # the same class b prior on y ~ 1 + x leaves the intercept alone
  f1 <- frm(bf(y ~ x), data = ri_data,
            prior = set_prior("normal(0, 0.1)", class = "b"))
  expect_gt(abs(fixef(f1)["Intercept", 1] - fe["Intercept", 1]),
            10 * fe["Intercept", 2])
  ps <- prior_summary(fit)
  expect_identical(ps$class, "b")
  tab <- as.data.frame(get_prior(bf(y ~ 0 + Intercept + x),
                                 data = ri_data))
  expect_false("Intercept" %in% tab$class)
  expect_true("Intercept" %in% tab$coef[tab$class == "b"])
  expect_error(
    frm(bf(y ~ 0 + Intercept + x), data = ri_data,
        prior = set_prior("normal(0, 1)", class = "Intercept")),
    "ordinary coefficient, as in brms: address it with class = \"b\"")
})

test_that("bf(center = FALSE) and lf(center = FALSE) are the same model", {
  p <- set_prior("normal(0, 0.1)", class = "b")
  a <- frm(bf(y ~ 0 + Intercept + x), data = ri_data, prior = p)
  b <- frm(bf(y ~ x, center = FALSE), data = ri_data, prior = p)
  expect_identical(fixef(a), fixef(b))
  # center = FALSE is the location formula's alone, as in brms
  tab <- as.data.frame(get_prior(bf(y ~ x, sigma ~ z, center = FALSE),
                                 data = ri_data))
  expect_setequal(tab$class[tab$dpar == "sigma" & tab$coef == ""],
                  c("Intercept", "b"))
  expect_false(any(tab$class == "Intercept" & tab$dpar == ""))
  tab <- as.data.frame(get_prior(bf(y ~ x) + lf(sigma ~ z, center = FALSE),
                                 data = ri_data))
  expect_false(any(tab$class == "Intercept" & tab$dpar == "sigma"))
  expect_true(any(tab$class == "Intercept" & tab$dpar == ""))
  expect_identical(bf(bf(y ~ x), center = FALSE)$center, FALSE)
  expect_error(bf(y ~ x, center = "no"), "`center` must be TRUE or FALSE")
  expect_error(lf(sigma ~ z, center = NA), "`center` must be TRUE or FALSE")
})

test_that("prediction, emmeans and the brms names need no Intercept column", {
  f0 <- frm(bf(y ~ 0 + Intercept + x + f), data = ri_data)
  f1 <- frm(bf(y ~ x + f), data = ri_data)
  nd <- data.frame(x = c(0, 1, 2), f = factor(c("a", "b", "d"),
                                              levels = letters[1:4]))
  expect_identical(fitted(f0, newdata = nd), fitted(f1, newdata = nd))
  expect_identical(variables(f0), variables(f1))
  expect_identical(hypothesis(f0, "Intercept = 0")$hypothesis,
                   hypothesis(f1, "Intercept = 0")$hypothesis)
  expect_identical(names(conditional_effects(f0)), c("x", "f"))
  skip_if_not_installed("emmeans")
  expect_equal(summary(emmeans::emmeans(f0, ~ f)),
               summary(emmeans::emmeans(f1, ~ f)))
})

test_that("coefficients come in brms's order", {
  # brms reorders its parameter BLOCKS, putting a centered intercept
  # first; a b_Intercept that is an element of the b vector stays where
  # its column is
  fit <- frm(bf(y ~ 0 + x + Intercept), data = ri_data)
  expect_identical(rownames(fixef(fit)), c("x", "Intercept"))
  fit <- frm(bf(y ~ x, sigma ~ 0 + Intercept + z), data = ri_data)
  expect_identical(rownames(fixef(fit)),
                   c("Intercept", "x", "sigma_Intercept", "sigma_z"))
  expect_identical(sigma(frm(bf(y ~ x, sigma ~ 0 + Intercept),
                             data = ri_data)),
                   sigma(frm(bf(y ~ x, sigma ~ 1), data = ri_data)))
})

test_that("the spellings brms refuses or does not reserve are refused", {
  expect_error(frm(bf(y ~ 1 + Intercept + x), data = ri_data),
               "reserved name of the intercept only in a population-level")
  expect_error(frm(bf(y ~ x + (0 + Intercept | f)), data = ri_data),
               "A group-level intercept is [(]1 [|] g[)]")
  expect_error(frm(bf(y ~ 0 + Intercept:x), data = ri_data),
               "Intercept:x is x")
  expect_error(frm(bf(y ~ 0 + Intercept + I(2 * Intercept)), data = ri_data),
               "only as a term of its own")
  expect_error(frm(bf(y ~ 0 + intercept + x), data = ri_data),
               "deprecated spelling of the reserved variable `Intercept`")
  d <- ri_data
  d$yo <- factor(cut(d$y, 4, labels = FALSE), ordered = TRUE)
  expect_error(frm(bf(yo ~ 0 + Intercept + x) + cumulative(), data = d),
               "whose thresholds take its place; brms refuses it too")
  d$Intercept <- 2
  expect_error(frm(bf(y ~ 0 + Intercept + x), data = d),
               "column `Intercept` that is not all ones")
  expect_error(frm(bf(y ~ x, sigma ~ 0 + Intercept), data = d),
               "column `Intercept` that is not all ones")
  # a column of ones is what brms would put there itself
  d$Intercept <- 1
  expect_identical(fixef(frm(bf(y ~ 0 + Intercept + x), data = d)),
                   fixef(frm(bf(y ~ 0 + Intercept + x), data = ri_data)))
})
