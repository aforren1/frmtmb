# The public accessors a custom density uses to read a distributional
# parameter on the linear-predictor scale.
#
# Two paths have to be covered by every one of these tests. ON the tape,
# build_objective() stores the linear predictor beside each dpar, and the
# accessor recovers the exact quantity from it. OFF the tape the entry is
# absent and the accessor falls back to the plain arithmetic. A family
# author who tests only off the tape never exercises the branch the
# accessors exist for, so both are asserted here side by side.

# The dpar list build_objective() hands a log-density, on the tape (with
# the linear predictor) and off it (without).
dp_on <- function(link, dpar, eta, extra = list()) {
  c(stats::setNames(list(link$linkinv(eta), eta),
                    c(dpar, paste0(".eta_", dpar))), extra)
}
dp_off <- function(link, dpar, eta, extra = list()) {
  c(stats::setNames(list(link$linkinv(eta)), dpar), extra)
}

lk_logit <- frmtmb:::get_link("logit")
lk_log <- frmtmb:::get_link("log")
lk_identity <- frmtmb:::get_link("identity")
lk_softplus <- frmtmb:::get_link("softplus")

test_that("the four accessors are exported", {
  # the point of the change: a family in another package can reach them
  for (nm in c("dpar_log", "dpar_log1m", "dpar_complement",
               "dpar_log_complement")) {
    expect_true(nm %in% getNamespaceExports("frmtmb"), info = nm)
    expect_true(is.function(get(nm, envir = asNamespace("frmtmb"))),
                info = nm)
  }
})

test_that("dpar_log is exact on the tape and plain off it", {
  for (eta in c(-700, -20, 0, 3, 40, 700)) {
    on <- frmtmb::dpar_log(dp_on(lk_log, "shape", eta), "shape", lk_log)
    expect_equal(on, eta, info = paste("eta", eta))
    off <- frmtmb::dpar_log(dp_off(lk_log, "shape", eta), "shape", lk_log)
    # off the tape the round trip is log(exp(eta)), which is the same
    # number wherever exp(eta) is representable and finite
    if (abs(eta) < 700) {
      expect_equal(off, eta, tolerance = 1e-12, info = paste("eta", eta))
    }
  }
  # exp(-700) is 9.86e-305 and its log is still -700; exp(-800) is 0 and
  # its log is -Inf. That is the gap the tape path closes.
  expect_true(is.infinite(frmtmb::dpar_log(
    dp_off(lk_log, "shape", -800), "shape", lk_log)))
  expect_equal(frmtmb::dpar_log(dp_on(lk_log, "shape", -800), "shape",
                                lk_log), -800)
})

test_that("dpar_log is exact for a unit-interval dpar too", {
  # a logit carries no log mean, only a log odds. Without that second
  # branch the same call is exact at eta = -800 on a log link and -Inf
  # on a logit, where plogis(-800) is exactly 0.
  expect_identical(log(stats::plogis(-800)), -Inf)
  # log(plogis(eta)) in the form that does not overflow at either end;
  # -log1p(exp(-eta)) is itself -Inf at eta = -800
  ref <- function(e) if (e > 0) -log1p(exp(-e)) else e - log1p(exp(e))
  for (eta in c(-800, -700, -40, 0, 40)) {
    v <- frmtmb::dpar_log(dp_on(lk_logit, "zi", eta), "zi", lk_logit)
    expect_equal(v, ref(eta), tolerance = 1e-12,
                 info = paste("eta", eta))
  }
})

test_that("dpar_log1m is the one-sided complement", {
  for (eta in c(-700, -40, 0, 40, 700)) {
    d <- dp_on(lk_logit, "coh", eta)
    v <- frmtmb::dpar_log1m(d, "coh", lk_logit)
    expect_true(is.finite(v), info = paste("eta", eta))
    expect_equal(v, -log1p(exp(eta)), tolerance = 1e-12,
                 info = paste("eta", eta))
    # and it agrees with the pair, which is the same log odds
    expect_identical(v, frmtmb::dpar_log_complement(d, "coh", lk_logit)$l1m,
                     info = paste("eta", eta))
  }
  # the subtraction it replaces is -Inf from eta = 36.7368005696771 up
  expect_true(is.infinite(log(1 - stats::plogis(40))))
  # off the tape it falls back, on log1p, and reaches the far end that
  # log(1 - p) cannot
  expect_identical(frmtmb::dpar_log1m(list(coh = 1e-17), "coh", lk_logit),
                   log1p(-1e-17))
  expect_identical(log(1 - 1e-17), 0)
})

test_that("dpar_log reads the dpar's OWN link, not an assumed log", {
  # a softplus predictor is not a log mean; reading it as one is a
  # different, wrong density
  eta <- 2
  v <- frmtmb::dpar_log(dp_on(lk_softplus, "shape", eta), "shape",
                        lk_softplus)
  expect_equal(v, log(log1p(exp(eta))))
  expect_false(isTRUE(all.equal(v, eta)))
})

test_that("dpar_log_complement stays finite where the subtraction dies", {
  # plogis(eta) is exactly 1 from about eta = 36.7368005696771 up, so
  # log(1 - plogis(eta)) is -Inf there and the density is unusable
  expect_true(is.infinite(log(1 - stats::plogis(40))))
  for (eta in c(-700, -40, -3, 0, 3, 40, 700)) {
    g <- frmtmb::dpar_log_complement(dp_on(lk_logit, "zi", eta), "zi",
                                     lk_logit)
    expect_true(is.finite(g$l), info = paste("l at eta", eta))
    expect_true(is.finite(g$l1m), info = paste("l1m at eta", eta))
    expect_equal(g$l1m, -log1p(exp(eta)), tolerance = 1e-12,
                 info = paste("eta", eta))
    expect_equal(g$l, -log1p(exp(-eta)), tolerance = 1e-12,
                 info = paste("eta", eta))
  }
  # off the tape the plain form is used, and is right where it is
  # representable
  g <- frmtmb::dpar_log_complement(dp_off(lk_logit, "zi", 1), "zi", lk_logit)
  p <- stats::plogis(1)
  expect_equal(g$l, log(p))
  expect_equal(g$l1m, log1p(-p))
  # and is NOT finite in the far tail, which is exactly the trap a family
  # author who tests only off the tape never sees
  expect_true(is.infinite(frmtmb::dpar_log_complement(
    dp_off(lk_logit, "zi", 40), "zi", lk_logit)$l1m))
})

test_that("dpar_log_complement reads the gate's OWN link", {
  # on an identity link there is no log odds to recover, so the pair has
  # to fall back rather than read a probability as a log odds
  d <- dp_on(lk_identity, "zi", 0.25)
  g <- frmtmb::dpar_log_complement(d, "zi", lk_identity)
  expect_equal(g$l, log(0.25))
  expect_equal(g$l1m, log(0.75))
})

test_that("dpar_complement is the exponential of the log pair", {
  for (eta in c(-40, 0, 40)) {
    d <- dp_on(lk_logit, "mu", eta)
    m <- frmtmb::dpar_complement(d, "mu", lk_logit)
    g <- frmtmb::dpar_log_complement(d, "mu", lk_logit)
    expect_equal(m$p, exp(g$l), info = paste("eta", eta))
    expect_equal(m$q, exp(g$l1m), info = paste("eta", eta))
  }
  # the complement is not zero where the subtraction makes it zero
  expect_identical(1 - stats::plogis(40), 0)
  expect_gt(frmtmb::dpar_complement(dp_on(lk_logit, "mu", 40), "mu",
                                    lk_logit)$q, 0)
  # off the tape it is the plain pair
  m <- frmtmb::dpar_complement(dp_off(lk_logit, "mu", 1), "mu", lk_logit)
  expect_equal(m$p, stats::plogis(1))
  expect_equal(m$q, 1 - stats::plogis(1))
})

test_that("a link may be named by a string", {
  # an extension declares links = list(coh = "logit") and never holds a
  # link object, so a name has to be accepted wherever an object is
  d <- dp_on(lk_logit, "coh", 40)
  expect_equal(frmtmb::dpar_log_complement(d, "coh", "logit"),
               frmtmb::dpar_log_complement(d, "coh", lk_logit))
  expect_equal(frmtmb::dpar_complement(d, "coh", "logit"),
               frmtmb::dpar_complement(d, "coh", lk_logit))
  expect_equal(frmtmb::dpar_log(dp_on(lk_log, "shape", 3), "shape", "log"),
               3)
  expect_equal(frmtmb::dpar_log1m(d, "coh", "logit"),
               frmtmb::dpar_log1m(d, "coh", lk_logit))
  expect_error(frmtmb::dpar_log(dp_on(lk_log, "shape", 3), "shape", "nope"),
               "Unknown link")
})

test_that("the accessors carry an exact gradient through the tape", {
  # value and gradient both, because the failure this closes is a NaN
  # gradient rather than only a NaN value
  f <- function(p) {
    g <- frmtmb::dpar_log_complement(
      dp_on(lk_logit, "zi", p[1]), "zi", lk_logit)
    -(g$l1m + 2 * g$l)
  }
  for (e0 in c(-40, 0, 40)) {
    tp <- RTMB::MakeTape(f, e0)
    expect_true(is.finite(tp(e0)), info = paste("value at", e0))
    gr <- as.numeric(tp$jacobian(e0))
    expect_true(is.finite(gr), info = paste("gradient at", e0))
    # d/deta [-(log(1-p) + 2 log p)] = p - 2(1 - p) = 3p - 2
    expect_equal(gr, 3 * stats::plogis(e0) - 2, tolerance = 1e-10,
                 info = paste("gradient at", e0))
  }
})

test_that("a custom family in a fit reaches the tape path", {
  # end to end: the accessors are useful only if build_objective()
  # really put the entries there, so this asserts through frm() rather
  # than through a hand-built dpar list
  seen <- new.env(parent = emptyenv())
  seen$eta <- FALSE
  fam <- frmtmb_family(
    "eta_probe", dpars = "mu", links = list(mu = "logit"),
    lpdf = function(y, dpars, aterms) {
      if (!is.null(dpars[[".eta_mu"]])) seen$eta <- TRUE
      g <- frmtmb::dpar_log_complement(dpars, "mu", "logit")
      y * g$l + (1 - y) * g$l1m
    },
    type = "discrete")
  set.seed(9)
  d <- data.frame(x = stats::rnorm(120))
  d$y <- stats::rbinom(120, 1, stats::plogis(0.5 + d$x))
  fit <- frm(bf(y ~ x), family = fam, data = d)
  expect_true(seen$eta)
  ref <- frm(bf(y ~ x), family = bernoulli(), data = d)
  expect_equal(as.numeric(stats::logLik(fit)),
               as.numeric(stats::logLik(ref)), tolerance = 1e-6)
  expect_equal(unlist(fixef(fit)$mu), unlist(fixef(ref)$mu),
               tolerance = 1e-5)
})

test_that("the fallback keeps log(1 - p) as p approaches zero", {
  # `1 - p` is exact for p >= 0.5, so the two spellings agree bit for
  # bit at the saturating end and the choice is free there.
  expect_identical(
    frmtmb::dpar_log_complement(list(zi = 0.999), "zi", lk_logit)$l1m,
    log(1 - 0.999))
  # They part company at the OTHER boundary, which a gate reaches just
  # as often: the subtraction rounds 1 - p to 1 and loses the value.
  expect_identical(log(1 - 1e-17), 0)
  g <- frmtmb::dpar_log_complement(list(zi = 1e-17), "zi", lk_logit)
  expect_identical(g$l1m, log1p(-1e-17))
  expect_lt(g$l1m, 0)
  # and the fallback still runs, and differentiates, ON the tape: a
  # link with no exact log-odds form (identity) reaches it there
  f <- function(x) {
    -frmtmb::dpar_log_complement(list(zi = x[1], .eta_zi = x[1]), "zi",
                                 lk_identity)$l1m
  }
  tp <- RTMB::MakeTape(f, 0.25)
  expect_equal(tp(0.25), -log1p(-0.25))
  expect_equal(as.numeric(tp$jacobian(0.25)), 1 / 0.75, tolerance = 1e-10)
})

test_that("a link object is held to the same contract as a link name", {
  # The trap this whole topic exists for, reached through the argument
  # rather than through the arithmetic: a list that is not an frmtmb
  # link has no `logit_eta`, so taking it on trust means falling back
  # and returning -Inf where the accessor exists to return -40. The
  # most likely such list is the one base R hands out, which spells the
  # derivative `mu.eta` where this contract spells it `mu_eta`.
  d40 <- dp_on(lk_logit, "zi", 40)
  expect_equal(frmtmb::dpar_log1m(d40, "zi", "logit"), -40,
               tolerance = 1e-12)
  expect_error(frmtmb::dpar_log1m(d40, "zi", stats::make.link("logit")),
               "mu_eta")
  expect_error(frmtmb::dpar_complement(d40, "zi", list(name = "logit")),
               "linkfun")
  # a genuine frmtmb link object still passes through
  expect_identical(frmtmb::dpar_log1m(d40, "zi", lk_logit),
                   frmtmb::dpar_log1m(d40, "zi", "logit"))
})

test_that("a dpar that is not in the list is named, not propagated", {
  # NULL propagates rather than stopping: `1 - NULL` is numeric(0), so
  # a misspelled dpar used to come back as an empty density term.
  d <- dp_on(lk_logit, "zi", 2)
  for (f in list(frmtmb::dpar_log, frmtmb::dpar_log1m,
                 frmtmb::dpar_log_complement, frmtmb::dpar_complement)) {
    expect_error(f(d, "z1", lk_logit), "z1")
  }
  # the reserved entries are not offered as candidates
  expect_error(frmtmb::dpar_complement(d, "z1", lk_logit), "carries: zi")
  # and the fallback path refuses in the same way, off the tape
  expect_error(frmtmb::dpar_complement(list(zi = 0.5), "z1", lk_logit),
               "z1")
})
