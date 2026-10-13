# Standard errors and warnings the SE check got wrong at 0.69.0, and
# what punch round 1 of lane seflags corrected in the lane's first
# rules: an estimable prediction flagged, a non-estimable hypothesis
# passed, an sd the likelihood ignores or reads only far away, curved
# ridges and a simplex face with silent SEs, and near-null directions
# (R/se-check.R). Lane seflags, dev/seflags-findings.md.

flags_conds <- function(expr) {
  w <- character()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}

flags_lost_phrase <- "Standard errors are not available"

# s(x)'s fixed column is exactly a + b x, so (Intercept), x and
# s(x).fx1 are aliased; the mean of x sets a, and with it the
# intercept's share of the flat direction (dev/seflags-alias.R)
flags_alias_data <- function(m, seed = 1) {
  set.seed(seed)
  n <- 200
  x <- stats::rnorm(n)
  x <- x - mean(x) + m
  d <- data.frame(x = x, g = factor(rep(1:10, n / 10)))
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + stats::rnorm(10, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.3)
  d
}

test_that("a prediction along an alias with a kept parameter is estimable", {
  # 0.69.0 review B1, built under the reference BLAS: the intercept's
  # share of the alias (about 1e-5) is below what the noise lets the
  # Hessian resolve, so 0.69.0 kept the intercept's SE; the null basis
  # the prediction test read dropped that share, the rest was not flat,
  # and all five predictions got NaN and a warning. The alias is the
  # design's and exact, so all three of its coefficients lose their SE
  # (punch round 1)
  d <- flags_alias_data(1e-5)
  fit <- suppressWarnings(frm(y ~ x + s(x) + (1 | g), data = d))
  lost <- sdr_of(fit)$se_lost
  expect_true(all(c("(Intercept)", "x", "s(x).fx1") %in% names(lost)))
  nd <- d[c(1, 50, 100, 150, 200), ]
  r <- flags_conds(fitted(fit, newdata = nd, re_formula = NA))
  expect_false(any(grepl("predictions move along", r$warnings,
                         fixed = TRUE)))
  se <- r$value[, "Est.Error"]
  expect_true(all(is.finite(se)))
  # the same model without the alias (s(x) holds the linear term) has
  # the same likelihood, so the same prediction standard errors
  ref <- suppressWarnings(frm(y ~ s(x) + (1 | g), data = d))
  se_ref <- fitted(ref, newdata = nd, re_formula = NA)[, "Est.Error"]
  expect_equal(unname(se), unname(se_ref), tolerance = 0.01)
  # the absent case: a quantity along the alias still has none
  r <- flags_conds(hypothesis(fit, "x = 0"))
  expect_true(is.nan(r$value$hypothesis$Est.Error))
  expect_true(any(grepl("hypotheses move along", r$warnings,
                        fixed = TRUE)))
})

# The alias of flags_alias_data() beside an identified near-collinear
# pair z2 = z1 + N(0, s) (the review's dev/seflags-rev-join.R)
flags_pair_data <- function(s, seed = 1) {
  set.seed(seed)
  n <- 200
  x <- stats::rnorm(n) + 0.1
  z1 <- stats::rnorm(n)
  d <- data.frame(x = x, z1 = z1, z2 = z1 + stats::rnorm(n, 0, s),
                  g = factor(rep(1:10, n / 10)))
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + 0.3 * d$z1 + 0.3 * d$z2 +
    stats::rnorm(10, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.3)
  d
}

test_that("a hypothesis off an exact alias has no SE, on it one", {
  # punch round 1, B1: with the pair beside the alias the noise bound on
  # the flat direction was 0.021, and b_x + t * bs_sx_1, estimable only
  # at t = b, got finite SEs and no warning up to 3 percent off b
  # (dev/seflags-b1ref.R). The design's alias is exact, so it is
  # tested exactly
  d <- flags_pair_data(1e-4)
  fit <- suppressWarnings(frm(y ~ x + s(x) + z1 + z2 + (1 | g), data = d))
  X <- fit$frame$linpreds[[1L]]$X
  b <- unname(stats::coef(stats::lm(as.numeric(X[, "s(x).fx1"]) ~
                                      d$x))[2L])
  hyp <- function(t) {
    flags_conds(hypothesis(fit, sprintf("b_x + %.15g * bs_sx_1 = 0", t),
                           class = NULL))
  }
  r <- hyp(b)
  expect_true(is.finite(r$value$hypothesis$Est.Error))
  expect_length(r$warnings, 0L)
  for (dl in c(1e-3, 1e-2)) {
    r <- hyp(b * (1 - dl))
    lab <- paste("delta", dl)
    expect_true(is.nan(r$value$hypothesis$Est.Error), info = lab)
    expect_true(any(grepl("hypotheses move along", r$warnings,
                          fixed = TRUE)), info = lab)
  }
})

test_that("an identified pair keeps its SEs beside an exact alias", {
  # punch round 1, B2: the rule that joined a confirmed direction to a
  # removed one it could not resolve took z1 and z2 as "flat" on 10 of
  # 15 such fits; they are identified (the fit without the alias, lm()
  # and lmer() report them), and the alias is still named
  d <- flags_pair_data(1e-5)
  fit <- suppressWarnings(frm(y ~ x + s(x) + z1 + z2 + (1 | g), data = d))
  ref <- suppressWarnings(frm(y ~ s(x) + z1 + z2 + (1 | g), data = d))
  lost <- sdr_of(fit)$se_lost
  expect_true(all(c("x", "s(x).fx1") %in% names(lost)))
  expect_false(any(c("z1", "z2") %in% names(lost)))
  se <- fixef(fit)[c("z1", "z2"), "Est.Error"]
  se_ref <- fixef(ref)[c("z1", "z2"), "Est.Error"]
  expect_equal(unname(se), unname(se_ref), tolerance = 1e-3)
  # the same alias across two nonlinear parameters, which no single
  # design shows, is named by the re-measured Hessian
  r <- flags_conds(frm(bf(y ~ a + b, a ~ 1 + x + z1 + z2 + (1 | g),
                          b ~ 1, nl = TRUE), data = flags_pair_data(3e-5, 2)))
  lost <- sdr_of(r$value)$se_lost
  expect_setequal(as.character(names(lost)),
                  c("a_(Intercept)", "b_(Intercept)"))
})

test_that("an sd read only far above its estimate is at its boundary", {
  # punch round 1, B3: y ~ a + 1e-4 * u with no group effect in the
  # data. The objective gains about 3e-7 moving the log sd down to -12
  # and changes by less than 1e-6 over the probe grid, but by 10 or more
  # at the log sd + 12; the first lane build called it "flat" or said
  # nothing (SE 778 to 995 on the log sd), 0.69.0 a boundary
  for (G in c(3, 10)) {
    set.seed(1)
    n <- 20 * G
    d <- data.frame(x = stats::rnorm(n), g = factor(rep(seq_len(G),
                                                        each = 20)))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(n, 0, 1)
    d$kk <- 1e-4
    r <- flags_conds(frm(bf(y ~ a + kk * u, a ~ 1 + x, u ~ 0 + (1 | g),
                            nl = TRUE), data = d))
    bnd <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"),
                  r$messages)
    expect_length(bnd, 1L)
    expect_identical(unname(sdr_of(r$value)$se_lost["theta_1"]),
                     "boundary", info = paste("G", G))
  }
})

flags_peak_data <- function(n_id = 10, seed = 5) {
  set.seed(seed)
  w <- seq(1, 45, by = 0.25)
  d <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    mu_i <- stats::rnorm(1, 10, 1.1)
    logf <- stats::rnorm(1, 1.2, 0.3) - 1.5 * log(w) +
      1.1 * exp(-0.5 * ((w - mu_i) / 1.6)^2)
    data.frame(id = i, w = w, logw = log(w),
               I = stats::rexp(length(w), rate = 1 / exp(logf)))
  }))
  d$id <- factor(d$id)
  d
}

test_that("an sd the likelihood does not read is flat, not at zero", {
  # 0.69.0 review, item 6: the bump is zero on every row, so the sd of
  # pk's group effect (0.22) does not enter the likelihood at any value;
  # 0.69.0 said "sd_id__pk_Intercept is at zero"
  d <- flags_peak_data()
  r <- flags_conds(frm(
    bf(I ~ apo - chi * logw + exp(lamp) * exp(-0.5 * ((w - pk) / exp(lsig))^2),
       apo ~ 1 + (1 | id), chi ~ 1, pk ~ 1 + (1 | id), lamp ~ 1, lsig ~ 1,
       nl = TRUE),
    family = exponential(link = "log"), data = d))
  fit <- r$value
  nm <- frmtmb:::outer_par_names(fit)
  j <- match("theta_2", nm)
  # the construction: the objective does not move with that sd
  p <- fit$opt$par
  f0 <- fit$obj$fn(p)
  for (h in c(-3, 1)) {
    q <- p
    q[j] <- q[j] + h
    expect_identical(fit$obj$fn(q), f0)
  }
  fit$obj$fn(p)
  bnd <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"), r$messages)
  expect_length(bnd, 0L)
  lost <- sdr_of(fit)$se_lost
  expect_identical(unname(lost["theta_2"]), "flat")
  w <- grep(flags_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "theta_2", fixed = TRUE)
})

flags_ridge_data <- function(seed) {
  set.seed(seed)
  n <- 60
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                   g = factor(rep(1:6, 10)))
  dd$y <- 3 + 0.5 * dd$x + 0.3 * dd$z + stats::rnorm(6, 0, 0.5)[dd$g] +
    stats::rnorm(n, 0, 0.4)
  dd
}

test_that("a curved ridge above the flat threshold loses its SEs", {
  # setier review m8: only a_(Intercept) + exp(b) is identified. The fit
  # stops within its tolerance of the ridge, where the Hessian along it
  # is 8.7e-7 (exact, seed 4) and 7.4e-7 (with (1 | g), seed 1), above
  # the flat threshold; both kept every standard error in silence
  # (dev/setier-rev-aexpb.R, dev/seflags-spectrum.R)
  dd <- flags_ridge_data(4)
  r <- flags_conds(frm(bf(y ~ a + exp(b), a ~ 1 + x, b ~ 1, nl = TRUE),
                       data = dd, start = list(beta = c(1, 0, -1))))
  # the construction: the likelihood does not move along the ridge
  fit <- r$value
  nm <- frmtmb:::outer_par_names(fit)
  ia <- match("a_(Intercept)", nm)
  ib <- match("b_(Intercept)", nm)
  p <- fit$opt$par
  q <- p
  q[ib] <- p[ib] + log(2)
  q[ia] <- p[ia] + exp(p[ib]) - exp(q[ib])
  f0 <- fit$obj$fn(p)
  expect_lt(abs(fit$obj$fn(q) - f0), 1e-9 * abs(f0))
  fit$obj$fn(p)
  lost <- sdr_of(fit)$se_lost
  expect_true(all(c("a_(Intercept)", "b_(Intercept)") %in% names(lost)))
  expect_false("a_x" %in% names(lost))
  expect_length(grep(flags_lost_phrase, r$warnings, fixed = TRUE), 1L)
  dd <- flags_ridge_data(1)
  r <- flags_conds(frm(bf(y ~ a + exp(b), a ~ 1 + x + (1 | g), b ~ 1,
                          nl = TRUE), data = dd,
                       start = list(beta = c(1, 0, -1))))
  lost <- sdr_of(r$value)$se_lost
  expect_true(all(c("a_(Intercept)", "b_(Intercept)") %in% names(lost)))
  se <- sqrt(diag(vcov(r$value, full = TRUE)))
  nm <- frmtmb:::outer_par_names(r$value)
  expect_true(is.nan(se[[match("b_(Intercept)", nm)]]))
})

test_that("a simplex face off the coordinates is found on every seed", {
  # optima's filing: with the third of four categories unobserved, two
  # steps of mo(inc) enter only as their sum, a flat direction that is
  # curved in the simplex chart; dev/optima-segap.R's seed 2 of
  # `ls ~ mo(inc)` and seed 3 of `ls ~ mo(inc) + age` kept every SE in
  # silence
  for (cs in list(list(ls ~ mo(inc), 2L), list(ls ~ mo(inc) + age, 3L))) {
    set.seed(cs[[2]])
    inc <- sample(c(0L, 1L, 3L), 100, TRUE)
    d <- data.frame(inc, ls = c(30, 60, 70, 75)[inc + 1L] +
                      stats::rnorm(100, sd = 7),
                    age = stats::rnorm(100, mean = 40, sd = 10))
    r <- flags_conds(frm(cs[[1]], data = d))
    lost <- sdr_of(r$value)$se_lost
    expect_setequal(as.character(names(lost)), c("zeta1_1", "zeta1_2"))
    expect_length(grep(flags_lost_phrase, r$warnings, fixed = TRUE), 1L)
  }
  # the absent case: every category observed keeps every SE in silence
  set.seed(2)
  inc <- sample(0:3, 100, TRUE)
  d <- data.frame(inc, ls = c(30, 60, 70, 75)[inc + 1L] +
                    stats::rnorm(100, sd = 7))
  r <- flags_conds(frm(ls ~ mo(inc), data = d))
  expect_length(sdr_of(r$value)$se_lost, 0L)
  expect_length(r$warnings, 0L)
})

test_that("a near-null direction next to a kept one keeps no shrunken SE", {
  skip_if_not_installed("lme4")
  # setier review m7: a degree-7 raw polynomial with (1 | g), seed 7. The
  # removed eigenvalue (8.9e-15 of the largest) sat 1.4e-11 from a kept
  # one; 0.69.0 named x3 to x5 alone and kept the others' standard
  # errors with neither near-null direction in them (the intercept's 330
  # against lmer()'s 2774), and the first lane build lost all eight as
  # "flat". The design is identified: re-measured by second differences
  # of the objective, both directions are real, and every coefficient's
  # SE is lmer()'s within 10 percent (dev/seflags-deg7diag.R; the review's
  # dev/seflags-rev-deg7fix.R gives 1.09 by the same measurement)
  set.seed(7)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  fit <- suppressWarnings(frm(fo, data = d))
  m <- allow_warnings(lme4::lmer(fo, data = d, REML = FALSE),
                      c("unidentifiable", "failed to converge"))
  ref <- sqrt(diag(as.matrix(stats::vcov(m))))
  se <- unname(fixef(fit)[, "Est.Error"])
  expect_length(sdr_of(fit)$se_lost, 0L)
  expect_true(all(is.finite(se)))
  expect_true(all(abs(se / ref - 1) < 0.15))
  vc <- sqrt(diag(vcov(fit, full = TRUE)))
  nm <- frmtmb:::outer_par_names(fit)
  expect_true(all(is.finite(vc[nm %in% c("sigma_(Intercept)", "theta_1")])))
})

test_that("a direction below double precision is ill-conditioned, not flat", {
  # A degree-8 raw polynomial with (1 | g), seed 3: its smallest
  # curvature is below machine epsilon of the largest. 0.69.0 called six
  # coefficients "flat" (the data do not determine them) and kept three.
  # The objective measures that curvature accurately, but one
  # double-precision covariance cannot hold it beside the others, and
  # punch round 2b, which kept it beside the matrix, gave emmeans,
  # marginaleffects and ordinal category SEs 0.3 to 255929 times the
  # true ones in silence (B6, B7). So the coefficients it carries have
  # no SE, nor does a prediction that reads them, and the reason says
  # the design is ill-conditioned, with the remedy, and that the fit
  # stopped short of the maximum (a Newton step gains 0.6 here)
  set.seed(3)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + (1 | g)
  # vcov() delivers a check that waited for the first SE use
  r <- flags_conds({
    m <- frm(fo, data = d)
    vcov(m)
    m
  })
  lost <- sdr_of(r$value)$se_lost
  expect_setequal(as.character(names(lost)),
                  c("(Intercept)", paste0("x", 1:8)))
  expect_true(all(lost == "illcond"))
  w <- grep(flags_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "the design is ill-conditioned", fixed = TRUE)
  expect_match(w, "poly()", fixed = TRUE)
  expect_match(w, "stopped short of the maximum", fixed = TRUE)
  expect_false(grepl("do not determine", w, fixed = TRUE))
  # B4: the first punch build gave the row at x = 2.5 an SE 0.31 times
  # the exact one, in silence; every row reads the lost coefficients
  nd <- data.frame(x = c(1, 1.5, 2, 2.5), g = d$g[1])
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  p <- flags_conds(fitted(r$value, newdata = nd, re_formula = NA))
  expect_true(all(is.nan(p$value[, "Est.Error"])))
  expect_true(any(grepl("predictions move along", p$warnings,
                        fixed = TRUE)))
})

test_that("terms beside an ill-conditioned polynomial keep their SEs", {
  # B5 of punch round 2: a variance-share rule took the SE of a
  # covariate z, a factor w and sigma beside a degree-8 raw polynomial,
  # which fixef() still showed. They do not read the lost coefficients,
  # so they keep their SEs: within 2 percent of the same likelihood
  # with poly(x, 8), well conditioned and at its maximum, or within 3
  # where the message says the fit stopped short of it (on a merge with
  # lane optrap's units these fits stop 1.1 to 1.5 short and sigma's SE
  # is 2.3 percent above the maximum's; dev/seflags-p2b-gain.R). Seeds 3
  # and 10 are the review's; under OpenBLAS's LAPACK this lane's
  # optimizer stops them at false convergence, which the SE check does
  # not read (filed), so there the test asks only that the fit said so.
  # Seed 4 converges on every build and asks both questions
  b5 <- function(seed) {
    set.seed(seed)
    d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                    z = stats::rnorm(200),
                    w = factor(sample(c("a", "b"), 200, TRUE)))
    d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
      stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
    for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
    raw <- flags_conds({
      m <- frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + w + (1 | g),
               data = d)
      vcov(m)
      m
    })
    ref <- frm(y ~ poly(x, 8) + z + w + (1 | g), data = d)
    conv <- raw$value$opt$convergence == 0
    expect_true(any(grepl(if (conv) "the design is ill-conditioned" else
      "did not report convergence", raw$warnings, fixed = TRUE)))
    if (conv) {
      lost <- sdr_of(raw$value)$se_lost
      expect_false(any(c("z", "wb") %in% names(lost)))
    }
    list(d = d, raw = raw$value, ref = ref, conv = conv,
         tol = if (any(grepl("stopped short", raw$warnings))) 0.03 else 0.02)
  }
  hyp <- function(f) {
    h <- flags_conds(hypothesis(f$raw, "b_wb = 0", class = NULL))
    expect_length(h$warnings, 0L)
    hr <- hypothesis(f$ref, "b_wb = 0", class = NULL)
    expect_lt(abs(h$value$hypothesis$Est.Error /
                    hr$hypothesis$Est.Error - 1), f$tol)
  }
  sig <- function(f) {
    nd <- data.frame(x = c(1, 1.5, 2), z = c(-1, 0, 1),
                     w = c("a", "b", "a"), g = f$d$g[1])
    for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
    s <- flags_conds(fitted(f$raw, newdata = nd, re_formula = NA,
                            dpar = "sigma"))
    expect_length(s$warnings, 0L)
    sr <- fitted(f$ref, newdata = nd, re_formula = NA, dpar = "sigma")
    expect_true(all(abs(s$value[, "Est.Error"] / sr[, "Est.Error"] - 1) <
                      f$tol))
  }
  f <- b5(3)
  if (f$conv) hyp(f)
  f <- b5(10)
  if (f$conv) sig(f)
  f <- b5(4)
  expect_true(f$conv)
  hyp(f)
  sig(f)
})

test_that("a prediction on a resolved near-null direction keeps its SE", {
  # the absent case of B4: degree 7 with (1 | g), seed 7, whose two
  # smallest directions the objective re-measures; fitted() rows in
  # and beyond the data's range are the exact prediction SE (whitened
  # QR at the fit's sigma and sd) within 1 percent. The first punch
  # build read solve(Q) when nothing was lost: 3.6 percent low at
  # x = 2.5, and 0 there on optrap's units
  set.seed(7)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  fit <- suppressWarnings(frm(fo, data = d))
  nd <- data.frame(x = c(1, 1.5, 2, 2.5), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  p <- flags_conds(fitted(fit, newdata = nd, re_formula = NA))
  expect_length(p$warnings, 0L)
  on <- frmtmb:::outer_par_names(fit)
  sg <- exp(fit$opt$par[grepl("^sigma", on)])
  tg <- exp(fit$opt$par[on == "theta_1"])
  X <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
  X0 <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, nd)
  Z <- stats::model.matrix(~ 0 + g, d)
  L <- chol(sg^2 * diag(nrow(X)) + tg^2 * Z %*% t(Z))
  q <- qr(backsolve(L, X, transpose = TRUE), tol = 1e-30)
  A <- backsolve(qr.R(q), t(X0[, q$pivot, drop = FALSE]), transpose = TRUE)
  ref <- sqrt(colSums(A^2))
  se <- unname(p$value[, "Est.Error"])
  expect_true(all(abs(se / ref - 1) < 0.01))
})
