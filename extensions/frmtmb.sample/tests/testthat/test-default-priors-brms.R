# frm_sample()'s default priors against brms 2.23.0's default_prior() on
# the same data (dev/correct-log/brms-priors.txt, script
# dev/correct-brms-probe.R). Through frmtmb.sample 0.9.0 a multivariate
# model got no defaults at all (rescor "(flat)" where brms has lkj(1)),
# the intercept location ignored an offset, and one class-wide sd default
# stood for brms's per-prefix sd rows.

dpb_rows <- function(tab) {
  tab <- as.data.frame(tab)
  stats::setNames(tab$prior, paste(tab$class, tab$coef, tab$group, tab$resp,
                                   tab$dpar, sep = "|"))
}

test_that("a multivariate model gets brms's rows, rescor lkj(1) included", {
  # brms's tests.priors.R "default_prior returns correct priors for
  # multivariate models", seed 1
  set.seed(1)
  dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), x = rnorm(10),
                    g = rep(1:2, 5))
  b15 <- bf(mvbind(y1, y2) ~ x + (x | ID1 | g)) + set_rescor(TRUE)
  r <- dpb_rows(default_prior(b15, dat, family = gaussian(),
                              route = "sample"))
  # brms 2.23.0's rows, as printed there
  want <- c("rescor||||" = "lkj(1)", "cor||||" = "lkj(1)",
            "Intercept|||y1|" = "student_t(3, 0.3, 2.5)",
            "Intercept|||y2|" = "student_t(3, 2, 2.5)",
            "sd|||y1|" = "student_t(3, 0, 2.5)",
            "sd|||y2|" = "student_t(3, 0, 2.5)",
            "sigma|||y1|" = "student_t(3, 0, 2.5)",
            "sigma|||y2|" = "student_t(3, 0, 2.5)",
            "b|||y1|" = "(flat)", "b|x||y1|" = "(flat)")
  expect_identical(r[names(want)], want)
})

test_that("the intercept location takes the mean offset off, as in brms", {
  # tests.priors.R "overall intercept priors are adjusted for the
  # intercept"
  dat <- data.frame(y = rep(c(1, 3), each = 5), off = 10)
  r <- dpb_rows(default_prior(y ~ 1 + offset(off), dat, route = "sample"))
  expect_identical(r[["Intercept||||"]], "student_t(3, -8, 2.5)")
  # a varying offset, printed to brms's digits
  set.seed(106)
  n <- 50
  d <- data.frame(x = rnorm(n), z = rnorm(n), o = rnorm(n, 3, 0.5),
                  o2 = rnorm(n, 0.4, 0.1))
  d$y <- 2 + d$o + d$x + rnorm(n)
  r <- dpb_rows(default_prior(y ~ x + offset(o), d, route = "sample"))
  expect_identical(r[["Intercept||||"]],
                   "student_t(3, 1.78834576782975, 2.5)")
  # an offset in sigma's formula does not move sigma's default
  r <- dpb_rows(default_prior(bf(y ~ x + offset(o), sigma ~ z + offset(o2)),
                              d, route = "sample"))
  expect_identical(r[["Intercept||||sigma"]], "student_t(3, 0, 2.5)")
  # log link: brms -0.900791310755016
  set.seed(107)
  dp <- data.frame(x = rnorm(n), expo = runif(n, 5, 50))
  dp$cnt <- rpois(n, dp$expo * exp(-1 + 0.3 * dp$x))
  r <- dpb_rows(default_prior(bf(cnt ~ x + offset(log(expo))) + poisson(),
                              dp, route = "sample"))
  expect_identical(r[["Intercept||||"]],
                   "student_t(3, -0.900791310755016, 2.5)")
  # and what is sampled is what is printed, to the 15 digits brms's
  # generated Stan code carries too
  uf <- frm(bf(y ~ x + offset(o)) + gaussian(), d, dry_run = "objective")
  s <- Filter(function(s) identical(s$class, "Intercept") &&
                !nzchar(s$dpar),
              unclass(frmtmb.sample:::default_priors_for(uf)))[[1L]]
  expect_identical(s$dist$location, as.numeric(as.character(
    round(stats::median(d$y), 1) - mean(d$o))))
})

dpb_beta_form <- function() {
  set.seed(104)
  n <- 120
  d <- data.frame(x = rnorm(n), g = factor(rep(1:8, length.out = n)))
  mu <- stats::plogis(0.2 + 0.4 * d$x + rnorm(8, 0, 0.4)[d$g])
  phi <- exp(2 + rnorm(8, 0, 0.3)[d$g])
  d$yb <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  list(d = d, form = bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta())
}

test_that("the sd default is written once per prefix", {
  fd <- dpb_beta_form()
  uf <- frm(fd$form, fd$d, dry_run = "objective")
  defs <- unclass(frmtmb.sample:::default_priors_for(uf))
  sdp <- Filter(function(s) identical(s$class, "sd"), defs)
  expect_setequal(vapply(sdp, `[[`, "", "dpar"), c("", "phi"))
  r <- dpb_rows(default_prior(fd$form, fd$d, route = "sample"))
  # brms 2.23.0: sd "" and sd phi both student_t(3, 0, 2.5)
  expect_identical(r[["sd||||"]], "student_t(3, 0, 2.5)")
  expect_identical(r[["sd||||phi"]], "student_t(3, 0, 2.5)")
})

# Its own block: on 0.61.0 the `sd phi` row above does not exist and the
# lookup ERRORS with "subscript out of bounds", so in one block these
# five assertions never ran on the base arm. They are about what a
# prior REACHES, which 0.61.0 can answer, so they are worth reaching.
test_that("a user sd prior reaches the location block only", {
  fd <- dpb_beta_form()
  uf <- frm(fd$form, fd$d, dry_run = "objective")
  # without a user sd prior every standard deviation keeps the density
  # the single class-wide default gave it before: nothing moves
  ri0 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf, NULL))$ri
  th <- Filter(function(z) identical(z$comp, "theta"), ri0$entries)
  expect_length(th, 2L)
  expect_true(all(vapply(th, function(z) z$dist$scale, 1) == 2.5))
  # a user sd prior with no dpar replaces the location parameter's
  # default and leaves phi's, as brms's stancode does
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(
    uf, set_prior("normal(0, 5)", class = "sd")))$ri
  th <- Filter(function(z) identical(z$comp, "theta"), ri$entries)
  kinds <- vapply(th, function(z) z$dist$kind, "")
  expect_setequal(kinds, c("normal", "t"))
})

test_that("prior_summary() of sampled multivariate draws lists brms's rows", {
  skip_on_cran()
  set.seed(105)
  n <- 80
  d <- data.frame(x = rnorm(n), g = factor(rep(1:8, length.out = n)))
  d$y1 <- 1 + d$x + rnorm(8)[d$g] + rnorm(n)
  d$y2 <- -1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(n)
  form <- bf(mvbind(y1, y2) ~ x + (1 | g)) + set_rescor(TRUE) + gaussian()
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, data = d, chains = 1, iter = 200, refresh = 0,
               seed = 5)))
  ps <- as.data.frame(prior_summary(ds))
  key <- paste(ps$class, ps$resp, sep = "|")
  got <- stats::setNames(ps$prior, key)
  # every row brms's default_prior() gives this model a density
  # (dev/correct-log/brms-priors.txt layout, this data's own numbers)
  tab <- as.data.frame(default_prior(form, d, route = "sample"))
  tab <- tab[tab$prior != "(flat)" & !nzchar(tab$coef) & !nzchar(tab$group),
             ]
  want <- stats::setNames(tab$prior, paste(tab$class, tab$resp, sep = "|"))
  expect_true(all(c("rescor|", "Intercept|y1", "Intercept|y2", "sd|y1",
                    "sd|y2", "sigma|y1", "sigma|y2") %in% names(want)))
  expect_identical(got[names(want)], want)
  expect_identical(got[["rescor|"]], "lkj(1)")
})
