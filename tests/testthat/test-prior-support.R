# A prior density with a bounded support keeps it on the tape, and a
# bounded density on an unbounded parameter warns as brms does
# (dev/reviews/2026-10-14-sampleprior.md, findings 2 and 5).

pst_tape <- function(dist) {
  RTMB::MakeTape(function(v) frmtmb:::prior_base_logdens(v, dist), 0)
}

test_that("each bounded density is NaN outside its support, on the tape", {
  # (shape - 1) * log(x) at shape 1 is a constant-zero operand the tape
  # folds (dev/rtmb-pitfalls.md, item 22), so through frmtmb 0.71.0
  # beta(1, 1) was 0 at -0.2 and 1.2 on the tape and gamma(1, 1) 0.5 at
  # -0.5; exponential() read no log at all and was finite below 0 on
  # doubles as well
  cases <- list(
    list(frmtmb:::parse_prior_dist("beta(1, 1)"), c(-0.2, 1.2), 0.3),
    list(frmtmb:::parse_prior_dist("beta(2, 1)"), c(-0.2, 1.2), 0.3),
    list(frmtmb:::parse_prior_dist("gamma(1, 1)"), -0.5, 0.7),
    list(frmtmb:::parse_prior_dist("gamma(2, 0.1)"), -0.5, 0.7),
    list(frmtmb:::parse_prior_dist("inv_gamma(0.4, 0.3)"), -0.5, 0.7),
    list(frmtmb:::parse_prior_dist("exponential(1)"), -0.5, 0.7),
    list(frmtmb:::parse_prior_dist("lognormal(0, 1)"), -0.5, 0.7),
    list(frmtmb:::parse_prior_dist("uniform(0, 2)"), c(-0.5, 2.5), 0.7))
  for (cs in cases) {
    lab <- frmtmb:::prior_dist_string(cs[[1]])
    tp <- pst_tape(cs[[1]])
    for (out in cs[[2]]) {
      # log() of a value outside warns on doubles; the tape is silent
      expect_true(is.nan(allow_warnings(
        frmtmb:::prior_base_logdens(out, cs[[1]]), "NaNs produced")),
                  label = paste(lab, "doubles at", out))
      expect_true(is.nan(tp(out)), label = paste(lab, "tape at", out))
    }
    expect_true(is.finite(tp(cs[[3]])), label = paste(lab, "inside"))
  }
})

test_that("the support is closed, as Stan's is", {
  # a value a link rounds onto the edge (plogis(40) is 1), or a mode on
  # it, keeps the density there; the round-1 guards made the edge NaN,
  # and an rdm() fit under ndt_group() whose ML ndt fraction is 1 then
  # could not start a chain (dev/sampleprior-rdm.R)
  edge <- list(
    list("uniform(0, 2)", c(0, 2), -log(2)),
    list("uniform(0.3, 0.7)", c(0.3, 0.7), -log(0.4)),
    list("beta(1, 1)", c(0, 1), 0),
    list("gamma(1, 1)", 0, 0),
    list("exponential(2)", 0, log(2)))
  for (cs in edge) {
    tp <- pst_tape(frmtmb:::parse_prior_dist(cs[[1]]))
    for (x in cs[[2]]) {
      expect_equal(tp(x), cs[[3]], tolerance = 1e-14,
                   label = paste(cs[[1]], "at", x))
      expect_true(all(is.finite(tp$jacobian(x))),
                  label = paste(cs[[1]], "gradient at", x))
    }
    # and one step past it is outside
    out <- c(min(cs[[2]]) - 1e-12, max(cs[[2]]) + 1e-12)
    if (!grepl("^(uniform|beta)", cs[[1]])) out <- out[1L]
    for (x in out) {
      expect_true(is.nan(tp(x)), label = paste(cs[[1]], "at", x))
    }
  }
})

test_that("inside the support the guarded densities are the densities", {
  # R's own densities, and the tape's gradient against a central
  # difference of them: the guard adds an exact 0 there
  ref <- list(
    "beta(1, 1)" = function(x) stats::dbeta(x, 1, 1, log = TRUE),
    "beta(2, 5)" = function(x) stats::dbeta(x, 2, 5, log = TRUE),
    "gamma(1, 1)" = function(x) stats::dgamma(x, 1, 1, log = TRUE),
    "gamma(0.01, 0.01)" = function(x) stats::dgamma(x, 0.01, 0.01,
                                                    log = TRUE),
    "exponential(2)" = function(x) stats::dexp(x, 2, log = TRUE),
    "lognormal(0.3, 1.2)" = function(x) stats::dlnorm(x, 0.3, 1.2,
                                                      log = TRUE),
    "inv_gamma(0.4, 0.3)" = function(x) {
      0.4 * log(0.3) - lgamma(0.4) - 1.4 * log(x) - 0.3 / x
    })
  for (nm in names(ref)) {
    dist <- frmtmb:::parse_prior_dist(nm)
    tp <- pst_tape(dist)
    for (x in c(0.05, 0.4, 0.93)) {
      expect_equal(frmtmb:::prior_base_logdens(x, dist), ref[[nm]](x),
                   tolerance = 1e-13, label = paste(nm, x))
      expect_equal(tp(x), frmtmb:::prior_base_logdens(x, dist),
                   tolerance = 1e-14)
      h <- 1e-6 * x
      fd <- (ref[[nm]](x + h) - ref[[nm]](x - h)) / (2 * h)
      g <- tp$jacobian(x)
      expect_equal(as.numeric(g), fd, tolerance = 1e-6,
                   label = paste(nm, "gradient", x))
    }
  }
})

test_that("a bounded density on an unbounded parameter warns as brms does", {
  # brms 2.23.0's check_prior_content(), measured with brms::stancode()
  # (dev/sampleprior-brmsunif.R): two warnings for uniform() on a
  # coefficient, one for gamma() there, the upper one only for uniform()
  # on sd or sigma, none for exponential() on sd or lognormal() on sigma,
  # and none once lb and ub are given. The check reads the prior a call
  # is given, as a whole (check_prior_slots(), which frm(), frm_sample(),
  # validate_prior() and update() run), not each set_prior() call
  warns <- function(pl) {
    w <- character(0)
    withCallingHandlers(check_prior_slots(pl), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    w
  }
  w <- warns(set_prior("uniform(2.2, 2.6)", class = "b", coef = "x"))
  expect_length(w, 2L)
  expect_match(w[1], "lower bounded prior on a parameter that has no",
               fixed = TRUE)
  expect_match(w[1], "b_x ~ uniform(2.2, 2.6)", fixed = TRUE)
  expect_match(w[2], "argument 'ub' of 'set_prior'", fixed = TRUE)
  expect_length(warns(set_prior("gamma(2, 1)", class = "b", coef = "x")), 1L)
  expect_length(warns(set_prior("beta(2, 2)", class = "Intercept")), 2L)
  w <- warns(set_prior("uniform(0, 10)", class = "sd"))
  expect_length(w, 1L)
  expect_match(w, "upper bounded", fixed = TRUE)
  expect_length(warns(set_prior("uniform(0, 10)", class = "sigma")), 1L)
  expect_length(warns(set_prior("exponential(1)", class = "sd")), 0L)
  expect_length(warns(set_prior("lognormal(0, 1)", class = "sigma")), 0L)
  expect_length(warns(set_prior("uniform(2.2, 2.6)", class = "b",
                                coef = "x", lb = 2.2, ub = 2.6)), 0L)
  # the absent case: a density with no bound, on any class
  expect_length(warns(set_prior("normal(0, 1)", class = "b")), 0L)
  # brms's bounds on a dpar's own class: beta() on zi, hu1 of a mixture
  expect_length(warns(set_prior("beta(2, 2)", class = "zi")), 0L)
  expect_length(warns(set_prior("beta(4, 16)", class = "hu1")), 0L)
})

test_that("a class row's bound covers its coefficients, as in brms", {
  # brms 2.23.0 reads the bounds of the combined prior, where a class
  # row's lb is the base bound of each coefficient under it, and does
  # not warn on this prior (dev/reviews/2026-10-14-sampleprior.md,
  # finding 10); frmtmb applies the class row's bound to each
  # coefficient too
  warns <- function(pl) {
    w <- character(0)
    withCallingHandlers(check_prior_slots(pl), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    w
  }
  expect_length(warns(set_prior("", class = "b", lb = 0) +
                        set_prior("gamma(2, 1)", class = "b", coef = "x")),
                0L)
  # the absent case: a class row with no bound does not cover them
  expect_length(warns(set_prior("normal(0, 5)", class = "b") +
                        set_prior("gamma(2, 1)", class = "b", coef = "x")),
                1L)
  # nor does a class row of another dpar
  expect_length(warns(set_prior("", class = "b", dpar = "sigma", lb = 0) +
                        set_prior("gamma(2, 1)", class = "b", coef = "x")),
                1L)
  # an lb on the class row leaves a two-sided density's ub to warn about
  w <- warns(set_prior("", class = "b", lb = 0) +
               set_prior("uniform(0, 3)", class = "b", coef = "x"))
  expect_length(w, 1L)
  expect_match(w, "upper bounded", fixed = TRUE)
  # and the warning comes from the call that uses the prior, once
  set.seed(20261024)
  d <- data.frame(x = stats::rnorm(30))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(30)
  w <- character(0)
  withCallingHandlers(
    frm(bf(y ~ x), data = d, dry_run = "objective",
        prior = set_prior("gamma(2, 1)", class = "b", coef = "x")),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_identical(sum(grepl("lower bounded prior", w, fixed = TRUE)), 1L)
})

test_that("the shared nu of a Student-t rescor model takes a prior", {
  # brms writes it with no resp; through frmtmb 0.71.0 it was refused
  # with resp (one nu across the responses) and without
  # (resp_missing_refusal())
  set.seed(20261114)
  n <- 60
  x <- stats::rnorm(n)
  d <- data.frame(x = x, y1 = 1 + x + stats::rt(n, 5),
                  y2 = 2 - x + stats::rt(n, 5))
  fit <- frm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = d,
             family = student(), dry_run = "objective")
  ri <- resolve_prior_input(fit, set_prior("gamma(2, 0.1)", class = "nu"))
  e <- Filter(function(z) identical(z$dist$kind, "gamma"), ri$entries)
  expect_length(e, 1L)
  expect_match(names(fit$frame$par_template[[e[[1]]$comp]])[e[[1]]$idx],
               "nu")
  tab <- validate_prior(set_prior("gamma(2, 0.1)", class = "nu"),
                        bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE),
                        data = d, family = student())
  expect_identical(tab$prior[tab$class == "nu"], "gamma(2, 0.1)")
  # with resp it stays refused, as brms requires
  expect_error(resolve_prior_input(fit, set_prior("gamma(2, 0.1)",
                                                  class = "nu",
                                                  resp = "y1")),
               "takes no resp")
})

test_that("a logit link's Jacobian is finite where plogis() rounds to 1", {
  # log(p (1 - p)) was -Inf past eta near 37, where 1 - p is 0, so a
  # natural-scale density on a dpar whose mode sits there (an rdm() ndt
  # fraction at 1, dev/sampleprior-rdm2.R) made the objective infinite
  lk <- frmtmb:::get_link("logit")
  for (eta in c(-50, -3, 0.4, 3, 40, 50)) {
    j <- frmtmb:::natural_log_jacobian(eta, lk)
    want <- stats::plogis(eta, log.p = TRUE) +
      stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)
    expect_equal(j, want, tolerance = 1e-13, label = paste("eta", eta))
  }
  tp <- RTMB::MakeTape(function(v) frmtmb:::natural_log_jacobian(v, lk), 0)
  expect_true(is.finite(tp(40)))
  expect_true(all(is.finite(tp$jacobian(40))))
})
