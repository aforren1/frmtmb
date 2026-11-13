# One rule for the user's own cold start (fit_assembled(), cold_start):
# the restart at a prior's location, the second reading of a gp() body and
# the mixture starts run on it alone, both candidates of the autoscale
# choice included; and a message that reports a point of the fit describes
# the fit frm() returns.

cs_conds <- function(expr, allowed = character()) {
  got <- list()
  val <- allow_warnings(withCallingHandlers(expr,
    warning = function(w) got[[length(got) + 1L]] <<- w,
    message = function(m) {
      got[[length(got) + 1L]] <<- m
      invokeRestart("muffleMessage")
    }), allowed = allowed)
  list(value = val, conds = got)
}
cs_has <- function(conds, cls) {
  any(vapply(conds, function(cnd) inherits(cnd, cls), NA))
}
cs_msgs <- function(conds) {
  # a message's text ends in a newline
  trimws(vapply(Filter(function(cnd) inherits(cnd, "message"), conds),
                conditionMessage, ""))
}
# the conflict of test-fit-restart.R: response near 1e6, slope 5e4, a
# normal(0, 1000) prior on the slope
cs_data <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(200))
  d$yl <- 1e6 * (1 + 0.05 * d$x) + rnorm(200, 0, 1e3)
  d
}
cs_prior <- function() {
  c(set_prior("normal(0, 1000)", class = "b"),
    set_prior("normal(1e6, 1e6)", class = "Intercept"))
}

test_that("the autoscale choice replays the messages of the fit it keeps", {
  # each candidate says something of its own point; the one kept is told
  # by the tier of mixture_starts(), which decides before the objective
  fake <- function(tag, tier, objective) {
    function() {
      frm_message("start of the ", tag, " fit",
                  class = "frmtmb_prior_multistart")
      list(opt = list(objective = objective,
                      mixture_starts = c(tier = tier)))
    }
  }
  run_with <- function(plain, seeded) {
    function(ctl, tpl) if (is.null(tpl)) plain() else seeded()
  }
  ctl <- frmtmb_control()
  got <- cs_conds(autoscale_choose(
    run_with(fake("plain", 2, 10), fake("seeded", 1, 20)), ctl,
    template = list()))
  expect_identical(got$value$opt$objective, 20)
  expect_identical(cs_msgs(got$conds), "start of the seeded fit")
  got <- cs_conds(autoscale_choose(
    run_with(fake("plain", 1, 10), fake("seeded", 2, 5)), ctl,
    template = list()))
  expect_identical(got$value$opt$objective, 10)
  expect_identical(cs_msgs(got$conds), "start of the plain fit")
})

test_that("both candidates of the autoscale choice take the cold start", {
  skip_on_cran()
  d <- cs_data(6201)
  calls <- 0L
  # one entry per mixture_starts() call that may run starts: its cold flag
  seen <- logical()
  orig_esc <- prior_conflict_escape
  orig_ms <- mixture_starts
  local_mocked_bindings(
    autoscale_prefit_verdict = function(...) "choose",
    prior_conflict_escape = function(...) {
      calls <<- calls + 1L
      orig_esc(...)
    },
    mixture_starts = function(fit, control, ...) {
      a <- list(...)
      if (!identical(control$starts, 0L)) seen <<- c(seen, isTRUE(a$cold))
      orig_ms(fit, control, ...)
    })
  got <- cs_conds(frm(yl ~ x, data = d, prior = cs_prior(),
                      control = frmtmb_control(autoscale = TRUE)),
                  allowed = ".")
  # the plain and the seeded candidate; the pre-fit runs neither search
  expect_identical(calls, 2L)
  expect_identical(seen, c(TRUE, TRUE))
  # the restart's message is said at most once, by the fit returned
  n <- sum(vapply(got$conds, inherits, NA, "frmtmb_prior_multistart"))
  expect_lte(n, 1L)
  expect_identical(n > 0L, !is.null(got$value$opt$prior_escape))
  # a refit is handed the fit's mode: no search runs
  calls <- 0L
  seen <- logical()
  rf <- cs_conds(refit(got$value, d$yl * (1 + 1e-12)), allowed = ".")
  expect_identical(calls, 0L)
  expect_false(any(seen))
})

test_that("a mixture start that wins takes the restart's message with it", {
  skip_on_cran()
  d <- cs_data(6201)
  # mixture_starts() hands back a start's optimum, which carries none of
  # the escape records of the optimum it replaced
  local_mocked_bindings(mixture_starts = function(fit, ...) {
    o <- fit$opt
    o$prior_escape <- NULL
    o$sd_escape <- NULL
    list(moved = TRUE, opt = o)
  })
  got <- cs_conds(frm(yl ~ x, data = d, prior = cs_prior()))
  expect_null(got$value$opt$prior_escape)
  expect_false(cs_has(got$conds, "frmtmb_prior_multistart"))
  # and with no start winning the restart reports its own point
  local_mocked_bindings(mixture_starts = function(...) NULL)
  got <- cs_conds(frm(yl ~ x, data = d, prior = cs_prior()))
  expect_false(is.null(got$value$opt$prior_escape))
  expect_true(cs_has(got$conds, "frmtmb_prior_multistart"))
})

# the data of the 0.82.0 release review's probe 2: y near 1e6 with a
# skewed residual, a group shift and a step in a monotonic predictor
cs_escape_data <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(200))
  e <- stats::rexp(200, 1 / 1e3) - 1e3
  d$y <- 1e6 * (1 + 0.05 * d$x) + e
  d$g <- factor(rep(1:20, each = 10))
  d$y2 <- d$y + stats::rnorm(20, 0, 30)[d$g]
  d$m <- sample(1:4, 200, TRUE)
  d$y3 <- d$y + 2e3 * c(0, 0.1, 0.9, 1)[d$m]
  d
}
cs_escape_prior <- function() {
  c(set_prior("normal(0, 1000)", class = "b"),
    set_prior("normal(1e6, 1e6)", class = "Intercept"))
}

test_that("an escape after the restart keeps the restart's message", {
  skip_on_cran()
  # each fit reaches the prior's mode by the restart and then a later
  # step climbs from it: mo_search() (mo(m), seed 6201), the stationary
  # escape (skew_normal(), seed 6206) and the sd escape ((1 | g), seed
  # 6206). That step's optimum carried no restart record, and the fit
  # returned the second mode in silence.
  d1 <- cs_escape_data(6201)
  d6 <- cs_escape_data(6206)
  fits <- list(
    mo = function() frm(y3 ~ x + mo(m), data = d1, prior = cs_escape_prior()),
    skew = function() frm(y ~ x, data = d6, family = skew_normal(),
                          prior = cs_escape_prior()),
    sd = function() frm(y2 ~ x + (1 | g), data = d6,
                        prior = cs_escape_prior()))
  won <- FALSE
  orig <- prior_conflict_escape
  local_mocked_bindings(prior_conflict_escape = function(...) {
    o <- orig(...)
    won <<- !is.null(o$prior_escape)
    o
  })
  for (nm in names(fits)) {
    won <- FALSE
    got <- cs_conds(fits[[nm]](), allowed = ".")
    # the message and the record go together whether the restart won or
    # not; on these data it wins
    expect_identical(cs_has(got$conds, "frmtmb_prior_multistart"), won,
                     info = nm)
    expect_identical(!is.null(got$value$opt$prior_escape), won, info = nm)
    expect_true(won, info = nm)
  }
})
