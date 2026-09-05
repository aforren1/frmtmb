# One simplex per mo() TERM, numbered as brms numbers its special terms.
#
# The frame used to key the simplex on the mo() VARIABLE, so
# `mo(x)` and `mo(x):z` in one formula shared a shape and the model had
# two fewer free parameters than brms's. The v0.18 NEWS called the
# sharing a brms convention; it was not, and the claim was withdrawn.

mo_terms_data <- function(seed = 21, n = 400) {
  set.seed(seed)
  inc <- sample(0:3, n, replace = TRUE)
  w <- sample(0:2, n, replace = TRUE)
  z <- rnorm(n)
  # the main effect rises early, the interaction late: two shapes a
  # single simplex cannot hold at once
  main <- c(0, 0.7, 0.9, 1)[inc + 1]
  slope <- c(0, 0.05, 0.15, 1)[inc + 1]
  data.frame(y = 1 + 2 * main + (0.3 + 1.5 * slope) * z + rnorm(n, 0, 0.6),
             inc = inc, w = w, z = z, z2 = rnorm(n))
}

zeta_names <- function(fit) {
  sort(grep("^zeta", names(fit$frame$par_template), value = TRUE))
}

test_that("an mo() variable in two terms gets two simplexes, one term one", {
  dd <- mo_terms_data()

  fit_int <- frm(bf(y ~ mo(inc) * z) + gaussian(), data = dd)
  expect_identical(zeta_names(fit_int), c("zeta1", "zeta2"))
  # inc runs 0..3, so D = 3 and the simplex is held as 2 free
  # softmax coordinates
  expect_length(fit_int$frame$par_template[["zeta1"]], 2L)
  expect_length(fit_int$frame$par_template[["zeta2"]], 2L)

  fit_add <- frm(bf(y ~ mo(inc) + z) + gaussian(), data = dd)
  expect_identical(zeta_names(fit_add), "zeta1")

  # the colon spelling is the same model as the star spelling minus the
  # main effect, so it too gets exactly one simplex
  fit_col <- frm(bf(y ~ mo(inc):z) + gaussian(), data = dd)
  expect_identical(zeta_names(fit_col), "zeta1")
})

test_that("splitting the simplex cannot lose likelihood, and gains it here", {
  dd <- mo_terms_data()
  fit <- frm(bf(y ~ mo(inc) * z) + gaussian(), data = dd)

  # the shared simplex the package used to build, fitted directly: the
  # same design with p[5:7] driving both monotonic columns. It is a
  # constrained submodel, so its optimum is a lower bound, and on data
  # whose two shapes differ it is a strictly worse one.
  step <- function(zr) {
    s <- exp(c(0, zr))
    c(0, cumsum(s / sum(s))) * 3
  }
  nll_shared <- function(p) {
    m <- step(p[5:7])[dd$inc + 1]
    -sum(stats::dnorm(dd$y,
                      p[1] + p[2] * dd$z + p[3] * m + p[4] * m * dd$z,
                      exp(p[8]), log = TRUE))
  }
  op <- stats::optim(c(1, 0.3, 0.5, 0.3, 0, 0, 0, log(0.6)), nll_shared,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_gte(as.numeric(logLik(fit)), -op$value)
  expect_gt(as.numeric(logLik(fit)), -op$value + 1)
  # and the extra shape is paid for in the degrees of freedom
  expect_identical(attr(logLik(fit), "df"),
                   attr(logLik(frm(bf(y ~ mo(inc) + z) + gaussian(),
                                   data = dd)), "df") + 3L)
})

test_that("simplexes are numbered in brms's special-term order", {
  skip_unless_brms()
  dd <- mo_terms_data()

  # brms enumerates through terms(), which lists every main effect
  # before any interaction whatever the written order. Our parser keeps
  # the written order and puts mo(x):z ahead of the mo(x) that mo(x) * z
  # implies, so this is the assertion that the frame re-sorts.
  check <- function(bform, fform) {
    sdat <- brms_standata(bform, data = dd, family = gaussian())
    fr <- frm(fform, data = dd, dry_run = "frame")
    mo <- fr$linpreds[["y.mu"]]$mo
    expect_identical(length(mo), as.integer(sdat$Imo))
    expect_identical(vapply(mo, function(m) m$D, integer(1)),
                     as.integer(sdat$Jmo))
    for (j in seq_along(mo)) {
      # Xmo_<j> is the category coding of the variable in brms's j-th
      # special term, which is what pins the ORDER rather than only the
      # count: mo(inc) and mo(w) have different codes
      expect_equal(mo[[j]]$codes, as.integer(sdat[[paste0("Xmo_", j)]]),
                   ignore_attr = TRUE)
      # and zeta<j> is brms's simo_<j>, holding its D - 1 free
      # coordinates against brms's D simplex entries
      expect_identical(mo[[j]]$zeta, paste0("zeta", j))
      expect_length(fr$par_template[[mo[[j]]$zeta]],
                    as.integer(sdat$Jmo[j]) - 1L)
    }
    vapply(mo, function(m) m$label, character(1))
  }

  expect_identical(check(brms::bf(y ~ mo(inc) * z),
                         bf(y ~ mo(inc) * z) + gaussian()),
                   c("moinc", "moinc:z"))
  expect_identical(check(brms::bf(y ~ z * mo(inc)),
                         bf(y ~ z * mo(inc)) + gaussian()),
                   c("moinc", "moinc:z"))
  # the written order puts the interaction first; brms still enumerates
  # the two main effects before it
  expect_identical(check(brms::bf(y ~ mo(inc) + mo(inc):z + mo(w)),
                         bf(y ~ mo(inc) + mo(inc):z + mo(w)) + gaussian()),
                   c("moinc", "mow", "moinc:z"))
  # different category counts, so Jmo alone identifies which term owns
  # which simplex
  expect_identical(check(brms::bf(y ~ mo(w) + mo(inc):z),
                         bf(y ~ mo(w) + mo(inc):z) + gaussian()),
                   c("mow", "moinc:z"))
  # TWO interactions. The sort leans on order() being stable to keep
  # the written order within one interaction order, and this is the
  # only shape where getting that tie wrong changes the answer: both
  # main effects first, then both interactions in the order written.
  expect_identical(check(brms::bf(y ~ mo(inc) * z + mo(w) * z2),
                         bf(y ~ mo(inc) * z + mo(w) * z2) + gaussian()),
                   c("moinc", "mow", "moinc:z", "mow:z2"))
})

test_that("brms names each simplex after its own term", {
  skip_unless_brms()
  dd <- mo_terms_data()

  # brms's simo coef is the term's b label with the index of the mo()
  # variable inside that term appended, and frmtmb allows one mo() per
  # term, so the index is always 1. Nothing in frmtmb emits a prior row
  # for a simplex, but the labels have to agree for one to be written.
  gp <- as.data.frame(brms::get_prior(brms::bf(y ~ mo(inc) * z),
                                      data = dd, family = gaussian()))
  expect_setequal(gp$coef[gp$class == "simo"],
                  c("moinc1", "moinc:z1"))

  fr <- frm(bf(y ~ mo(inc) * z) + gaussian(), data = dd,
            dry_run = "frame")
  mo <- fr$linpreds[["y.mu"]]$mo
  expect_identical(paste0(vapply(mo, function(m) m$label, character(1)), 1L),
                   c("moinc1", "moinc:z1"))
})

test_that("mo() and mi() two-way interactions still fit together", {
  set.seed(31)
  n <- 250
  dd <- data.frame(inc = sample(0:3, n, replace = TRUE), z = rnorm(n))
  dd$x <- rnorm(n, 0.5 + 0.8 * dd$z, 0.7)
  dd$y <- 1 + c(0, 0.6, 0.9, 1)[dd$inc + 1] * 2 + 0.5 * dd$x +
    0.3 * dd$z + rnorm(n, 0, 0.7)
  dd$x[sample(n, 40)] <- NA

  fit <- frm(bf(y ~ mo(inc) * z + mi(x) * z) + bf(x | mi() ~ z) +
               gaussian(), data = dd)
  # the mo() terms take two simplexes; the mi() terms take none
  expect_identical(zeta_names(fit), c("zeta1", "zeta2"))
  # mi() terms keep the written order, which puts the interaction
  # first. They carry no simplex, so no zeta is numbered against brms.
  # Their POSITION is not free of brms, though: brms puts mo() and mi()
  # terms in one bsp, interleaved in terms() order, so it wants
  # mo(inc) | mi(x) | mo(inc):z | mi(x):z where this reports
  # moinc | moinc:z | mix:z | mix. Nothing asserts that today (no tier
  # row uses mi()), which is why the mi() list is left unsorted.
  expect_setequal(names(fixef(fit)$y_mu),
                  c("(Intercept)", "z", "moinc", "moinc:z", "mix",
                    "mix:z"))
  expect_true(is.finite(as.numeric(logLik(fit))))
  nd <- data.frame(inc = 0:3, z = 0, x = 0)
  expect_true(all(is.finite(predict(fit, newdata = nd, resp = "y"))))
})

test_that("the tier translates an ordinal fit, whose zetas start at 2", {
  skip_unless_brms()
  dd <- mo_terms_data()
  dd$yo <- as.integer(cut(dd$y, quantile(dd$y, 0:3 / 3),
                          include.lowest = TRUE))

  # the family claims the first extras slot (tau_raw), so the monotonic
  # simplexes are zeta2 and zeta3 while brms still calls them simo_1
  # and simo_2. The translator has to map by POSITION in the mo list;
  # reading the trailing integer would look for a zeta1 that is the
  # thresholds, or is absent.
  fit <- frm(bf(yo ~ mo(inc) * z) + cumulative(), data = dd)
  expect_identical(grep("^zeta", names(fit$estimates), value = TRUE),
                   c("zeta2", "zeta3"))

  bform <- brms::bf(yo ~ mo(inc) * z)
  prior <- brms_flat_prior(bform, data = dd, family = brms::cumulative())
  code <- brms::make_stancode(bform, data = dd,
                              family = brms::cumulative(), prior = prior)
  sdat <- brms_standata(bform, data = dd, family = brms::cumulative(),
                        prior = prior)
  pars <- stan_pars_from_fit(fit, sdat, code)

  expect_true(all(c("simo_1", "simo_2") %in% names(pars)))
  simplex <- function(z) {
    v <- exp(c(0, fit$estimates[[z]]))
    v / sum(v)
  }
  expect_equal(pars[["simo_1"]], simplex("zeta2"))
  expect_equal(pars[["simo_2"]], simplex("zeta3"))
  expect_length(pars[["simo_1"]], sdat$Jmo[[1]])
})

test_that("the tier refuses a monotonic simplex outside mu", {
  skip_unless_brms()
  dd <- mo_terms_data()

  # brms spells this simo_sigma_1. Its trailing integer would alias
  # onto mu's sequence and hand over mu's simplex, at mu's length, with
  # no error; no row fits mo() outside mu, so the rule refuses by name
  # rather than guess.
  fit <- frm(bf(y ~ mo(inc), sigma ~ mo(w)) + gaussian(), data = dd)
  bform <- brms::bf(y ~ mo(inc), sigma ~ mo(w))
  prior <- brms_flat_prior(bform, data = dd, family = gaussian())
  code <- brms::make_stancode(bform, data = dd, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = dd, family = gaussian(),
                        prior = prior)
  expect_error(stan_pars_from_fit(fit, sdat, code),
               "monotonic simplex outside")
})
