# Which standard deviations a class "sd" prior reaches, against brms
# 2.23.0's stancode on the same models (dev/correct-log/brms-priors.txt
# and brms-priors2.txt, script dev/correct-brms-probe.R). brms keys its
# sd rows by response, distributional parameter and nonlinear
# parameter, and an empty field there means the empty one, not "every":
#
#   bf(yb ~ x + (1 | g), phi ~ (1 | g)), Beta(), sd group = "g":
#     lprior += normal_lpdf(sd_1 | 0, 5)            the mu block only
#     lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)    phi keeps its default
#   mvbind(y1, y2) ~ x + (1 | g), sd group = "g", no resp: refused
#   nonlinear a ~ 1 + (1 | g), sd with no nlpar: refused
#
# Through 0.61.0 the first reached the phi block too and the second and
# third were accepted and broadcast to every block.

sds_reached <- function(fit, pl) {
  e <- frmtmb:::resolve_prior_input(fit, pl)$entries
  as.numeric(sort(unlist(lapply(e, function(z) {
    if (identical(z$comp, "theta")) z$idx
  }))))
}

sds_beta_data <- function() {
  set.seed(104)
  n <- 120
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:8, length.out = n)))
  mu <- stats::plogis(0.2 + 0.4 * d$x + stats::rnorm(8, 0, 0.4)[d$g])
  phi <- exp(2 + stats::rnorm(8, 0, 0.3)[d$g])
  d$yb <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  d
}

test_that("a class sd prior with no dpar stays on the location blocks", {
  d <- sds_beta_data()
  fit <- suppressWarnings(frm(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
                              data = d))
  # which theta is which block: the mu block's sd, then phi's
  th <- vapply(fit$frame$re_blocks, function(b) as.numeric(b$theta_idx[1]), 1)
  dp <- vapply(fit$frame$re_blocks, function(b) b$dpar, "")
  mu_th <- th[dp == "mu"]
  phi_th <- th[dp == "phi"]
  expect_length(mu_th, 1L)
  expect_length(phi_th, 1L)
  expect_identical(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                              group = "g")), mu_th)
  expect_identical(sds_reached(fit, set_prior("normal(0, 5)", class = "sd")),
                   mu_th)
  expect_identical(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                              dpar = "phi")), phi_th)
  expect_identical(sds_reached(fit, set_prior("normal(0, 5)", class = "sd") +
                                 set_prior("normal(0, 2)", class = "sd",
                                           dpar = "phi")),
                   sort(c(mu_th, phi_th)))
  # the table says the same thing: a row per prefix, the phi one default
  vp <- validate_prior(set_prior("normal(0, 5)", class = "sd", group = "g"),
                       bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
                       data = d)
  sd <- as.data.frame(vp)[vp$class == "sd", ]
  expect_identical(sd$prior[sd$group == "g" & sd$dpar == ""],
                   "normal(0, 5)")
  expect_identical(sd$prior[sd$group == "g" & sd$dpar == "phi"], "(flat)")
  expect_identical(sort(sd$dpar[sd$group == ""]), c("", "phi"))
})

test_that("a block across mu and sigma follows brms's per-prefix rule", {
  # brms: (1 | q | g) in mu and sigma is ONE block, sd_1[1] mu and
  # sd_1[2] sigma. sd group g: normal on the whole sd_1. sd dpar sigma:
  # normal on sd_1[2] only. So on a spanning block an empty field
  # reaches both, and a more specific specification takes its own.
  set.seed(101)
  n <- 120
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:10, length.out = n)))
  d$y <- 1 + d$x + stats::rnorm(10)[d$g] +
    stats::rnorm(n, 0, exp(0.2 + stats::rnorm(10, 0, 0.3)[d$g]))
  fit <- frm(bf(y ~ x + (1 | q | g), sigma ~ (1 | q | g)) + gaussian(),
             data = d)
  expect_length(fit$frame$re_blocks, 1L)
  bk <- fit$frame$re_blocks[[1]]
  sd_th <- as.numeric(bk$theta_idx[frmtmb:::block_sd_idx(bk)])
  expect_length(sd_th, 2L)
  expect_identical(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                              group = "g")), sd_th)
  expect_identical(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                              dpar = "sigma")), sd_th[2])
  e <- frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 5)", class = "sd") +
      set_prior("normal(0, 2)", class = "sd", dpar = "sigma"))$entries
  sc <- vapply(e, function(z) z$dist$scale %||% NA_real_, 1)
  ix <- vapply(e, function(z) as.numeric(z$idx[1]), 1)
  expect_identical(sc[ix == sd_th[1]], 5)
  expect_identical(sc[ix == sd_th[2]], 2)
})

test_that("a multivariate sd prior needs its response, as in brms", {
  set.seed(105)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:8, length.out = n)))
  d$y1 <- 1 + d$x + stats::rnorm(8)[d$g] + stats::rnorm(n)
  d$y2 <- -1 + 0.5 * d$x + stats::rnorm(8)[d$g] + stats::rnorm(n)
  fit <- frm(bf(mvbind(y1, y2) ~ x + (1 | g)) + set_rescor(TRUE) +
               gaussian(), data = d)
  expect_error(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                          group = "g")),
               "No random-effect SDs match", class = "frmtmb_error")
  expect_error(sds_reached(fit, set_prior("normal(0, 5)", class = "sd")),
               "resp = \"y1\"", class = "frmtmb_error")
  expect_error(frm(bf(mvbind(y1, y2) ~ x + (1 | g)) + set_rescor(TRUE) +
                     gaussian(), data = d,
                   prior = set_prior("normal(0, 5)", class = "sd")),
               "No random-effect SDs match", class = "frmtmb_error")
  r1 <- sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                   group = "g", resp = "y1"))
  r2 <- sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                   group = "g", resp = "y2"))
  expect_length(r1, 1L)
  expect_length(r2, 1L)
  expect_false(identical(r1, r2))
})

test_that("a nonlinear parameter's sd prior needs its nlpar, as in brms", {
  set.seed(111)
  n <- 120
  d <- data.frame(x = stats::runif(n, 0, 3),
                  g = factor(rep(1:8, length.out = n)))
  d$y <- (2 + stats::rnorm(8, 0, 0.3)[d$g]) * exp(-0.5 * d$x) +
    stats::rnorm(n, 0, 0.1)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
               gaussian(), data = d, dry_run = "objective")
  expect_error(sds_reached(fit, set_prior("normal(0, 5)", class = "sd")),
               "nlpar = \"a\"", class = "frmtmb_error")
  expect_error(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                          group = "g")),
               "No random-effect SDs match", class = "frmtmb_error")
  expect_length(sds_reached(fit, set_prior("normal(0, 5)", class = "sd",
                                           nlpar = "a")), 1L)
})
