# frm_ode() against the core's nonlinear-formula surface. This lived in
# frmtmb's own test-nlf.R until the solver moved out; it stays with the
# solver, because what it checks is that frm_ode() composes, not that
# nlf() works.

test_that("frm_ode() composes inside an nlf() body", {
  skip_if_not_installed("RTMBode")
  pk_dyn <- function(t, y, p) {
    list(c(-p[1] * y[1], p[1] * y[1] - p[2] / p[3] * y[2]))
  }
  set.seed(2026)
  tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  id <- factor(rep(seq_len(4), each = length(tt)))
  time <- rep(tt, 4)
  ka <- exp(stats::rnorm(4, 0, 0.3))[as.integer(id)]
  ke <- exp(log(0.2) + stats::rnorm(4, 0, 0.25))[as.integer(id)]
  mu <- 100 * ka / (10 * (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  d <- data.frame(id = id, time = time, dose = 100,
                  conc = mu + stats::rnorm(length(mu), 0, 0.3))
  st <- list(beta = c(0, log(0.25), log(8)))

  direct <- frm(bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0),
                                  times = time,
                                  parms = list(exp(lka), exp(lke),
                                               exp(lV)),
                                  group = id, output = 2L),
                   lka ~ 1, lke ~ 1, lV ~ 1, nl = TRUE) + gaussian(),
                data = d, start = st)
  vianlf <- frm(bf(conc ~ pk, nl = TRUE) +
                  nlf(pk ~ frm_ode(pk_dyn, init = list(dose, 0),
                                   times = time,
                                   parms = list(exp(lka), exp(lke),
                                                exp(lV)),
                                   group = id, output = 2L)) +
                  lf(lka ~ 1, lke ~ 1, lV ~ 1) + gaussian(),
                data = d, start = st)
  expect_equal(as.numeric(logLik(vianlf)), as.numeric(logLik(direct)),
               tolerance = 1e-10)
  expect_equal(predict(vianlf), predict(direct), tolerance = 1e-10)

  # the within-group constancy check reads every body, not just mu's
  d$phase <- factor(ifelse(d$time <= 2, "early", "late"))
  expect_error(
    frm(bf(conc ~ pk, nl = TRUE) +
          nlf(pk ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                           parms = list(exp(lka), exp(lke), exp(lV)),
                           group = id, output = 2L)) +
          lf(lka ~ 1, lke ~ 1 + phase, lV ~ 1) + gaussian(),
        data = d, dry_run = "frame",
        start = list(beta = c(0, log(0.25), 0, log(8)))),
    "not constant within 'id'"
  )
})

test_that("the frame check is registered with the core", {
  # the seam itself: if .onLoad() did not register, the constancy
  # refusal above would pass silently and the fit would be wrong
  expect_true(any(vapply(
    frmtmb:::frmtmb_frame_checks$fns,
    function(f) identical(environmentName(environment(f)), "frmtmb.ode"),
    logical(1)
  )))
})
