# A nonlinear body is a scope of its own: the functions it CALLS come
# from RTMB before anywhere else.
#
# Before this, a bare pnorm() in an `nl` body resolved lexically to
# stats::pnorm(), which broke the tape or - with RTMB attached by
# accident - quietly worked, so the same formula meant different things
# in different sessions. The body now carries RTMB's replacements for
# the base and stats functions it calls (`nl_rtmb_shadow`, R/ad-env.R).
#
# This file pins four things: the shadow set is the mechanical
# collision set and every member of it is numerically transparent; the
# bare spelling and the `RTMB::` spelling are the same fit bit for bit;
# RTMB wins over the formula environment while everything else in that
# environment still resolves; and the taped and numeric evaluations of
# one body agree, which is what the transparency audit is for.

prep_dd <- function(n = 200, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, -2, 2))
  d$y <- 1.4 * stats::pnorm(0.8 * d$x - 0.3) + stats::rnorm(n, 0, 0.15)
  d
}

probit_bf <- function(fn = quote(pnorm)) {
  # the same body with the head swapped, so the two spellings differ in
  # nothing else
  body <- bquote(a * .(fn)(b * x + c0))
  form <- stats::as.formula(call("~", quote(y), body))
  bf(form, a ~ 1, b ~ 1, c0 ~ 1, nl = TRUE)
}

# --- the shadow set is derived, not asserted by taste ----------------

test_that("the shadow set is the RTMB/base/stats collision, less qchisq", {
  collide <- sort(intersect(getNamespaceExports("RTMB"),
                            c(ls(baseenv()),
                              getNamespaceExports("stats"))))
  shipped <- frmtmb:::nl_rtmb_shadow
  # a name RTMB stopped exporting would be shadowed to nothing, and a
  # new colliding export would silently stay lexical: both are changes
  # to the promise, so both fail here rather than in a user's formula
  expect_setequal(shipped, setdiff(collide, "qchisq"))
  expect_false("qchisq" %in% shipped)
  expect_true(all(c("pnorm", "qnorm", "dnorm", "dgamma", "pgamma",
                    "qgamma", "dbeta", "pbeta", "dbinom", "dpois",
                    "dnbinom", "plogis", "qlogis", "besselK") %in%
                    shipped))
})

test_that("every shadowed name is numerically transparent", {
  # the numeric paths - simulate(), predict(), the plain-numeric
  # objective - re-run the same body off the tape. A shadowed name that
  # returned a different number there would split the two paths without
  # saying so, so transparency is the whole safety argument.
  q  <- c(-2.5, -0.7, 0, 0.3, 1.25, 3.1)
  p  <- c(0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999)
  xp <- c(0.05, 0.4, 1, 2.3, 7.7)
  u  <- c(0.02, 0.3, 0.6, 0.97)
  k  <- c(0, 1, 2, 5, 9)
  M  <- matrix(c(2, 0.4, 0.1, 0.4, 3, 0.2, 0.1, 0.2, 1.5), 3, 3)
  V  <- c(1.3, -0.6, 2.2)

  probes <- list(
    dnorm = function(f) list(f(q, 0.4, 1.7), f(q, 0.4, 1.7, log = TRUE)),
    pnorm = function(f) list(f(q, 0.4, 1.7), f(q, 0.4, 1.7, log.p = TRUE),
                             f(q, 0.4, 1.7, lower.tail = FALSE)),
    qnorm = function(f) list(f(p, 0.4, 1.7), f(p, 0.4, 1.7,
                                               lower.tail = FALSE)),
    dgamma = function(f) list(f(xp, 2.2, 1.4), f(xp, 2.2, 1.4, log = TRUE)),
    pgamma = function(f) list(f(xp, 2.2, 1.4), f(xp, 2.2, 1.4,
                                                 lower.tail = FALSE)),
    qgamma = function(f) list(f(p, 2.2, 1.4)),
    dbeta = function(f) list(f(u, 2, 3), f(u, 2, 3, log = TRUE)),
    pbeta = function(f) list(f(u, 2, 3)),
    qbeta = function(f) list(f(p, 2, 3)),
    dbinom = function(f) list(f(k, 9, 0.3), f(k, 9, 0.3, log = TRUE)),
    pbinom = function(f) list(f(k, 9, 0.3)),
    dpois = function(f) list(f(k, 2.4), f(k, 2.4, log = TRUE)),
    ppois = function(f) list(f(k, 2.4)),
    dnbinom = function(f) list(f(k, 3, 0.4), f(k, 3, mu = 2.5, log = TRUE)),
    pnbinom = function(f) list(f(k, 3, 0.4)),
    dexp = function(f) list(f(xp, 1.6), f(xp, 1.6, log = TRUE)),
    pexp = function(f) list(f(xp, 1.6)),
    qexp = function(f) list(f(p, 1.6)),
    dlnorm = function(f) list(f(xp, 0.2, 1.1), f(xp, 0.2, 1.1, log = TRUE)),
    dlogis = function(f) list(f(q, 0.3, 1.2), f(q, 0.3, 1.2, log = TRUE)),
    plogis = function(f) list(f(q, 0.3, 1.2), f(q, 0.3, 1.2, log.p = TRUE)),
    qlogis = function(f) list(f(p, 0.3, 1.2)),
    dcauchy = function(f) list(f(q, 0.3, 1.2), f(q, 0.3, 1.2, log = TRUE)),
    dchisq = function(f) list(f(xp, 4), f(xp, 4, log = TRUE)),
    pchisq = function(f) list(f(xp, 4)),
    dt = function(f) list(f(q, 7), f(q, 7, log = TRUE)),
    df = function(f) list(f(xp, 4, 9), f(xp, 4, 9, log = TRUE)),
    dweibull = function(f) list(f(xp, 1.8, 2.2), f(xp, 1.8, 2.2,
                                                   log = TRUE)),
    pweibull = function(f) list(f(xp, 1.8, 2.2)),
    qweibull = function(f) list(f(p, 1.8, 2.2)),
    dmultinom = function(f) list(f(c(2, 3, 5), prob = c(0.2, 0.3, 0.5)),
                                 f(c(2, 3, 5), prob = c(0.2, 0.3, 0.5),
                                   log = TRUE)),
    besselK = function(f) list(f(xp, 0.5), f(xp, 1.5), f(xp, 0.5, TRUE)),
    besselI = function(f) list(f(xp, 0.5), f(xp, 1.5)),
    besselJ = function(f) list(f(xp, 0.5), f(xp, 1.5)),
    besselY = function(f) list(f(xp, 0.5), f(xp, 1.5)),
    lbeta = function(f) list(f(xp, 2.2), f(2.2, xp)),
    cov2cor = function(f) list(f(M)),
    colSums = function(f) list(f(M)),
    rowSums = function(f) list(f(M)),
    diag = function(f) list(f(M), f(3), f(V)),
    matrix = function(f) list(f(1:6, 2, 3), f(1:6, 2, 3, byrow = TRUE)),
    solve = function(f) list(f(M), f(M, V)),
    eigen = function(f) list(f(M, symmetric = TRUE)$values),
    svd = function(f) list(f(M)$d),
    fft = function(f) list(f(V)),
    order = function(f) list(f(c(3, 1, 2)), f(c(3, 1, 2),
                                              decreasing = TRUE)),
    sort = function(f) list(f(c(3, 1, 2)), f(c(3, 1, 2),
                                             decreasing = TRUE)),
    ifelse = function(f) list(f(q > 0, q, -q)),
    findInterval = function(f) list(f(q, c(-1, 0, 1))),
    apply = function(f) list(f(M, 1, sum), f(M, 2, sum)),
    sapply = function(f) list(f(1:4, function(i) i^2)),
    Vectorize = function(f) list(f(function(a, b) a + b)(1:3, 4:6)),
    integrate = function(f) list(f(function(z) z^2, 0, 1)$value),
    uniroot = function(f) list(f(function(z) z^2 - 2, c(0, 4))$root),
    splinefun = function(f) list(f(1:5, c(1, 4, 9, 16, 25))(2.5))
  )

  # every shipped name is probed: an unprobed one would be a promise
  # made without evidence
  expect_setequal(names(probes), frmtmb:::nl_rtmb_shadow)

  st_ns <- asNamespace("stats")
  for (nm in frmtmb:::nl_rtmb_shadow) {
    orig <- if (exists(nm, envir = st_ns, inherits = FALSE)) {
      get(nm, envir = st_ns)
    } else {
      get(nm, envir = baseenv())
    }
    # tolerance = 0 : the same function, not merely a close one
    expect_equal(probes[[nm]](frmtmb:::nl_shadow_fun(nm)),
                 probes[[nm]](orig), tolerance = 0,
                 info = paste("shadowed name:", nm))
  }
})

test_that("qchisq is excluded because RTMB reimplements it", {
  # kept as the reason, not as taste: if RTMB ever makes qchisq exact,
  # this fails and the name can join the set
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_false(isTRUE(all.equal(RTMB::qchisq(p, 4), stats::qchisq(p, 4),
                                tolerance = 0)))
  expect_equal(RTMB::qchisq(p, 4), stats::qchisq(p, 4), tolerance = 1e-12)
})

# --- called names are shadowed, read names are not -------------------

test_that("nl_shadow_for shadows what a body calls, not what it reads", {
  expect_identical(frmtmb:::nl_shadow_for(quote(a * pnorm(b * x))),
                   "pnorm")
  expect_setequal(frmtmb:::nl_shadow_for(quote(qgamma(pnorm(z), s, r))),
                  c("qgamma", "pnorm"))
  # a READ name keeps its meaning: `df` is far more often a data frame
  # than the F density, and shadowing it would break a body that works
  expect_identical(frmtmb:::nl_shadow_for(quote(b0 * nrow(df))),
                   character(0))
  expect_identical(frmtmb:::nl_shadow_for(quote(b0 * order)), character(0))
  # both ways in one body is an ambiguity no formula should contain, so
  # the conservative half wins
  expect_identical(frmtmb:::nl_shadow_for(quote(df(x, 4, 9) + df$k)),
                   character(0))
  # a name that is not RTMB's is never touched
  expect_identical(frmtmb:::nl_shadow_for(quote(my_helper(x))),
                   character(0))
  expect_identical(frmtmb:::nl_shadow_for(NULL), character(0))
})

test_that("the enclosure shadows only the called names", {
  e <- frmtmb:::ad_overload_env(baseenv(), quote(pnorm(x) + my_helper(x)))
  expect_setequal(ls(e, all.names = TRUE),
                  c("c", "[<-", "diag<-", "pnorm"))
  expect_identical(get("pnorm", envir = e),
                   getExportedValue("RTMB", "pnorm"))
  # the formula environment is the parent, never replaced
  expect_identical(parent.env(e), baseenv())
  # and no body means no shadowing, which is what every non-nl caller
  # of this enclosure gets
  expect_setequal(ls(frmtmb:::ad_overload_env(baseenv()), all.names = TRUE),
                  c("c", "[<-", "diag<-"))
})

# --- the two spellings are one fit -----------------------------------

test_that("a bare pnorm body is the RTMB:: body bit for bit", {
  d <- prep_dd()
  st <- list(beta = c(1, 1, 0))
  bare <- frm(probit_bf(quote(pnorm)), family = gaussian(), data = d,
              start = st)
  qual <- frm(probit_bf(quote(RTMB::pnorm)), family = gaussian(),
              data = d, start = st)

  expect_identical(as.numeric(logLik(bare)), as.numeric(logLik(qual)))
  expect_identical(coef(bare), coef(qual))
  expect_identical(unname(fitted(bare)), unname(fitted(qual)))
  # and it is a real fit, not two matching failures
  expect_equal(unname(fixef(bare)$a[["(Intercept)"]]), 1.4,
               tolerance = 0.1)
  expect_equal(unname(fixef(bare)$b[["(Intercept)"]]), 0.8,
               tolerance = 0.1)
})

test_that("qgamma and plogis bodies tape bare, and match the prefix", {
  set.seed(11)
  n <- 150
  d <- data.frame(x = stats::runif(n, -1.5, 1.5))
  d$y <- 2 * stats::plogis(1.3 * d$x) + stats::rnorm(n, 0, 0.1)
  st <- list(beta = c(1, 1))
  b1 <- frm(bf(y ~ a * plogis(b * x), a ~ 1, b ~ 1, nl = TRUE),
            family = gaussian(), data = d, start = st)
  q1 <- frm(bf(y ~ a * RTMB::plogis(b * x), a ~ 1, b ~ 1, nl = TRUE),
            family = gaussian(), data = d, start = st)
  expect_identical(as.numeric(logLik(b1)), as.numeric(logLik(q1)))

  # qgamma is the trap the audit recorded: it takes a RATE third, and
  # the point here is only that the bare spelling reaches it
  set.seed(12)
  d2 <- data.frame(z = stats::rnorm(n))
  d2$y <- stats::qgamma(stats::pnorm(d2$z), 3, 2) + stats::rnorm(n, 0, 0.05)
  st2 <- list(beta = c(3, 2))
  b2 <- frm(bf(y ~ qgamma(pnorm(z), sh, rt), sh ~ 1, rt ~ 1, nl = TRUE),
            family = gaussian(), data = d2, start = st2)
  q2 <- frm(bf(y ~ RTMB::qgamma(RTMB::pnorm(z), sh, rt), sh ~ 1, rt ~ 1,
               nl = TRUE),
            family = gaussian(), data = d2, start = st2)
  expect_identical(as.numeric(logLik(b2)), as.numeric(logLik(q2)))
  expect_equal(unname(fixef(b2)$sh[["(Intercept)"]]), 3, tolerance = 0.15)
})

test_that("an nlf() body is shadowed the same way as an inline one", {
  d <- prep_dd()
  st <- list(beta = c(1, 1, 0))
  inline <- frm(probit_bf(quote(pnorm)), family = gaussian(), data = d,
                start = st)
  composed <- frm(bf(y ~ m, nl = TRUE) +
                    nlf(m ~ a * pnorm(b * x + c0)) +
                    lf(a ~ 1, b ~ 1, c0 ~ 1),
                  family = gaussian(), data = d, start = st)
  expect_identical(as.numeric(logLik(inline)),
                   as.numeric(logLik(composed)))
})

# --- the gradient, not just the value --------------------------------

test_that("a bare-pnorm body has the gradient finite differences give", {
  d <- prep_dd()
  fit <- frm(probit_bf(quote(pnorm)), family = gaussian(), data = d,
             start = list(beta = c(1, 1, 0)))
  # away from the optimum, where a wrong derivative cannot hide
  pv <- fit$opt$par + c(0.05, -0.07, 0.03, 0.02)
  an <- as.numeric(fit$obj$gr(pv))
  fd <- vapply(seq_along(pv), function(i) {
    h <- 1e-5 * max(1, abs(pv[i]))
    up <- pv; up[i] <- up[i] + h
    dn <- pv; dn[i] <- dn[i] - h
    (fit$obj$fn(up) - fit$obj$fn(dn)) / (2 * h)
  }, numeric(1))
  expect_equal(an, fd, tolerance = 1e-5)
})

# --- the taped and the numeric evaluation are one function -----------

test_that("a bare-pnorm body gives one logLik taped and numeric", {
  d <- prep_dd()
  fit <- frm(probit_bf(quote(pnorm)), family = gaussian(), data = d,
             start = list(beta = c(1, 1, 0)))
  # the plain-numeric closure: no MakeADFun, so the body runs off the
  # tape through the same enclosure
  numeric_nll <- frmtmb:::build_objective(fit$frame)(fit$estimates)
  expect_equal(as.numeric(numeric_nll), -as.numeric(logLik(fit)),
               tolerance = 1e-10)
  # and the taped objective agrees with the reported logLik
  expect_equal(fit$obj$fn(fit$opt$par), -as.numeric(logLik(fit)),
               tolerance = 1e-9, ignore_attr = TRUE)
  # eval_dpars() and the newdata route are the other two numeric paths
  expect_equal(unname(fitted(fit)),
               unname(predict(fit, type = "response")), tolerance = 1e-12)
  expect_equal(unname(predict(fit, newdata = d)), unname(predict(fit)),
               tolerance = 1e-12)
})

# --- precedence: RTMB wins, everything else still resolves -----------

test_that("a user-defined pnorm is NOT picked up inside a body", {
  d <- prep_dd()
  # the body is a language with its own vocabulary, so this one loses.
  # It would flatten the model if it were found, which is why the fit
  # below is the same fit as with no local pnorm at all.
  pnorm <- function(q, mean = 0, sd = 1, ...) rep(0.5, length(q))
  st <- list(beta = c(1, 1, 0))
  shadowed <- frm(probit_bf(quote(pnorm)), family = gaussian(), data = d,
                  start = st)
  reference <- frm(probit_bf(quote(RTMB::pnorm)), family = gaussian(),
                   data = d, start = st)
  expect_identical(as.numeric(logLik(shadowed)),
                   as.numeric(logLik(reference)))
  # the local definition is real, and still means itself outside a body
  expect_equal(pnorm(c(-1, 0, 1)), c(0.5, 0.5, 0.5))
})

test_that("stats:: is the escape hatch a body can still reach", {
  d <- prep_dd()
  st <- list(beta = c(1, 1, 0))
  # the strongest proof that a qualified name is honored rather than
  # rewritten: stats::pnorm cannot tape an advector, so asking for it
  # fails exactly as it always did. The shadowing reaches bare names
  # only.
  expect_error(frm(probit_bf(quote(stats::pnorm)), family = gaussian(),
                   data = d, start = st),
               "Non-numeric argument to mathematical function")

  # and where stats:: is usable it is used: qnorm(0.75) is a constant,
  # so the fit is the bare-pnorm fit with `a` absorbing the factor
  scaled <- frm(bf(y ~ a * pnorm(b * x + c0) * stats::qnorm(0.75),
                   a ~ 1, b ~ 1, c0 ~ 1, nl = TRUE),
                family = gaussian(), data = d, start = st)
  plain <- frm(probit_bf(quote(pnorm)), family = gaussian(), data = d,
               start = st)
  expect_equal(as.numeric(logLik(scaled)), as.numeric(logLik(plain)),
               tolerance = 1e-8)
  expect_equal(unname(fixef(scaled)$a[["(Intercept)"]]) *
                 stats::qnorm(0.75),
               unname(fixef(plain)$a[["(Intercept)"]]), tolerance = 1e-5)
})

test_that("a user helper in the formula environment is still found", {
  # the interposition must shadow names, not replace the enclosure
  d <- prep_dd()
  squash <- function(u) u / (1 + abs(u))
  fit <- frm(bf(y ~ a * squash(b * x), a ~ 1, b ~ 1, nl = TRUE),
             family = gaussian(), data = d, start = list(beta = c(1, 1)))
  expect_s3_class(fit, "frmtmb_fit")
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("a read name keeps the user's object even when RTMB has it", {
  # `df` is in the shadow set, and this body READS it. Binding the F
  # density here would turn a working fit into "non-numeric argument".
  d <- prep_dd()
  df <- data.frame(k = c(2, 3, 5))
  scale_by <- function(tab, v) v * nrow(tab)
  fit <- frm(bf(y ~ b0 * scale_by(df, x), b0 ~ 1, nl = TRUE),
             family = gaussian(), data = d, start = list(beta = 1))
  expect_s3_class(fit, "frmtmb_fit")
  # the same body with the object written inline is the same fit
  inl <- frm(bf(y ~ b0 * scale_by(data.frame(k = c(2, 3, 5)), x), b0 ~ 1,
                nl = TRUE),
             family = gaussian(), data = d, start = list(beta = 1))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(inl)),
               tolerance = 1e-12)
})

# --- the acceptance case ---------------------------------------------

test_that("a response-preparation model fits with no RTMB:: anywhere", {
  # Two Gaussian preparation processes race a deadline: the habitual
  # response is ready first and is usually wrong, the goal-directed one
  # is ready later and is usually right, and before either is ready the
  # answer is a guess. That gives the dip in accuracy at intermediate
  # preparation times which the two CDFs and the mixture arithmetic
  # below reproduce. Identity link, because the body already returns a
  # probability.
  set.seed(7)
  n <- 3000
  pt <- stats::runif(n, 0, 0.6)
  m1 <- 0.20; m2 <- 0.32; sdev <- 0.045
  f1 <- stats::pnorm(pt, m1, sdev)
  f2 <- stats::pnorm(pt, m2, sdev)
  p <- f2 * 0.95 + f1 * (1 - f2) * 0.2 + (1 - f1) * (1 - f2) * 0.5
  d <- data.frame(pt = pt, y = stats::rbinom(n, 1, p))

  # every start comes from a prior location, and no name is qualified
  pr <- prior(normal(0.15, 1), nlpar = "m1") +
    prior(normal(0.40, 1), nlpar = "m2") +
    prior(normal(-3, 1), nlpar = "ls")
  form <- bf(
    y ~ pnorm(pt, m2, exp(ls)) * 0.95 +
      pnorm(pt, m1, exp(ls)) * (1 - pnorm(pt, m2, exp(ls))) * 0.2 +
      (1 - pnorm(pt, m1, exp(ls))) * (1 - pnorm(pt, m2, exp(ls))) * 0.5,
    m1 ~ 1, m2 ~ 1, ls ~ 1, nl = TRUE)

  fit <- suppressMessages(
    frm(form, family = bernoulli(link = "identity"), data = d, prior = pr))

  expect_s3_class(fit, "frmtmb_fit")
  # absolute bounds, because these are times in seconds and a relative
  # tolerance on a 45 ms standard deviation says nothing useful
  expect_lt(abs(unname(fixef(fit)$m1[["(Intercept)"]]) - m1), 0.025)
  expect_lt(abs(unname(fixef(fit)$m2[["(Intercept)"]]) - m2), 0.015)
  expect_lt(abs(exp(unname(fixef(fit)$ls[["(Intercept)"]])) - sdev), 0.006)

  # the fitted probability is a probability on every row, and the two
  # numeric routes through the body agree with the taped one. The
  # likelihood-only closure is checked in the gaussian test above
  # instead: a prior is a penalty on this objective, so here
  # build_objective(frame) and logLik() differ by exactly that penalty.
  fv <- fitted(fit)
  expect_true(all(fv > 0 & fv < 1))
  expect_equal(fit$obj$fn(fit$opt$par), -as.numeric(logLik(fit)),
               tolerance = 1e-9, ignore_attr = TRUE)
  expect_equal(unname(fv), unname(predict(fit, type = "response")),
               tolerance = 1e-12)
  expect_equal(unname(predict(fit, newdata = d)), unname(predict(fit)),
               tolerance = 1e-12)
})
