# The CDF slot of frmtmb_family() is `cdf`. Its former name, `lcdf`, read
# as a log CDF beside `lccdf` and `cdf_log` while it returned
# probabilities, so it is refused on every route a family reaches core by,
# with an error that names `cdf`. No alias: a family that filled the old
# name would otherwise lose its CDF without a word.
#
# dev/cdfname-rename.R skips this file: it must spell the old name.

pois_family <- function(...) {
  frmtmb_family(
    "pois_cdfname", dpars = "mu", links = list(mu = "log"),
    type = "discrete", accepts_aterms = c("weights", "cens", "trunc"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars[["mu"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y)), ...)
}
pois_cdf <- function(q, dpars, aterms) RTMB::ppois(q, dpars[["mu"]])

# a family as an extension that still fills the old name would hand it
# over: built, then edited
old_named <- function() {
  f <- pois_family()
  f[["lcdf"]] <- pois_cdf
  f
}

cens_data <- function() {
  set.seed(3088)
  d <- data.frame(y = stats::rpois(80, 4))
  d$cc <- ifelse(d$y > 6, "right", "none")
  d$y <- pmin(d$y, 6)
  d
}

renamed <- function(expr) {
  expect_error(expr, "`cdf`", class = "frmtmb_cdf_slot_renamed")
}

test_that("the slot is `cdf` and returns probabilities", {
  f <- pois_family(cdf = pois_cdf)
  expect_true(is.function(f$cdf))
  expect_false("lcdf" %in% names(f))
  expect_equal(f$cdf(3, list(mu = 4), list()), stats::ppois(3, 4))
  # a core family fills the same slot
  expect_equal(as_frmtmb_family(poisson())$cdf(3, list(mu = 4), list()),
               stats::ppois(3, 4))
})

test_that("frmtmb_family(lcdf = ) is refused, naming `cdf`", {
  renamed(pois_family(lcdf = pois_cdf))
  # any other stray argument is an unused argument, not the rename
  err <- tryCatch(pois_family(cdff = pois_cdf), error = identity)
  expect_false(inherits(err, "frmtmb_cdf_slot_renamed"))
  expect_match(conditionMessage(err), "unused argument: `cdff`",
               fixed = TRUE)
})

test_that("a family that stores `lcdf` is refused by frm() and bf()", {
  d <- cens_data()
  renamed(frm(y | cens(cc) ~ 1, family = old_named(), data = d))
  renamed(frm(y ~ 1, family = old_named(), data = d))
  renamed(bf(y | cens(cc) ~ 1, family = old_named()))
  renamed(bf(y | cens(cc) ~ 1) + old_named())
  renamed(as_frmtmb_family(old_named()))
  renamed(as_frmtmb_family(old_named))
})

test_that("a mixture component that stores `lcdf` is refused", {
  renamed(mixture(old_named(), poisson()))
})

test_that("frame assembly refuses one that skipped as_frmtmb_family()", {
  d <- cens_data()
  b <- bf(y | cens(cc) ~ 1) + pois_family(cdf = pois_cdf)
  b$family[["lcdf"]] <- pois_cdf
  renamed(frm(b, data = d))
})

test_that("a family_finalize() that adds `lcdf` is refused", {
  d <- cens_data()
  f <- pois_family(cdf = pois_cdf, family_finalize = function(fam, y, av) {
    fam[["lcdf"]] <- fam[["cdf"]]
    fam
  })
  renamed(frm(y ~ 1, family = f, data = d))
})

test_that("reading `$lcdf` is an error naming `cdf`, not NULL", {
  renamed(as_frmtmb_family(poisson())$lcdf)
  # a name that is no field at all is still NULL
  expect_null(as_frmtmb_family(poisson())$no_such_slot)
})

test_that("an absent `cdf` element reads as NULL, not as `cdf_log`", {
  # `cdf` is a prefix of `cdf_log`, and `$` refuses a prefix; assigning
  # NULL drops the element, as an extension's seam may
  f <- pois_family()
  f[["cdf"]] <- NULL
  expect_false("cdf" %in% names(f))
  expect_null(f$cdf)
  expect_null(f$lccdf)
  # a name frmtmb_family() does not define is still refused as a prefix
  expect_error(f$cdf_l, "partial-matched `cdf_log`", fixed = TRUE)
})

test_that("unused unnamed arguments are counted", {
  err <- tryCatch(frmtmb:::refuse_family_dots(1, 2, x = 3),
                  error = identity)
  expect_match(conditionMessage(err), "unused arguments: `x`, 2 unnamed",
               fixed = TRUE)
})

# A fit saved before the rename, read back: every family object in it
# stores the slot under the old name
old_shape <- function(x) {
  if (inherits(x, "frmtmb_family")) {
    names(x)[names(x) == "cdf"] <- "lcdf"
    return(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    for (i in seq_along(x)) {
      if (is.list(x[[i]])) x[i] <- list(old_shape(x[[i]]))
    }
  }
  x
}

test_that("an old-shaped fit is refused where a method reads the CDF", {
  set.seed(3088)
  d <- data.frame(x = stats::rnorm(120))
  d$y <- stats::rweibull(120, shape = 1.5, scale = exp(0.3 + 0.4 * d$x))
  fit <- frm(y ~ x, family = weibull(), data = d)
  # the guard's condition absent: the same fit as built reads its CDF
  r <- residuals(fit, type = "osa")[, "Estimate"]
  expect_true(all(is.finite(r)))
  old <- old_shape(fit)
  # PIN: the method fell back to the one-step integral of the density,
  # with integrate() errors and a NaN, rather than naming the slot
  renamed(residuals(old, type = "osa"))
})

test_that("an old-shaped truncated fit is refused where draws invert", {
  # a truncated draw reads the CDF to invert it
  set.seed(3089)
  dt <- data.frame(y = pmin(stats::rpois(80, 4), 8))
  ft <- old_shape(frm(y | trunc(ub = 8) ~ 1, family = poisson(), data = dt))
  renamed(simulate(ft, nsim = 1, seed = 1))
})

test_that("with the slot named `cdf`, censoring fits as before", {
  # the guard's condition absent: the same family, slot spelled right,
  # matches core's poisson() on the same censored rows
  d <- cens_data()
  fit <- frm(y | cens(cc) ~ 1, family = pois_family(cdf = pois_cdf),
             data = d)
  ref <- frm(y | cens(cc) ~ 1, family = poisson(), data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-10)
})
