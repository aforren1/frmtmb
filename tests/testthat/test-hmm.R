## hmm(): first-class hidden Markov models.
##
## The reference implementations are depmixS4 (gaussian, poisson and
## categorical emissions, with and without transition covariates),
## hmmTMB (the stationary initial distribution and random effects), and
## a plain numeric forward / forward-backward / Viterbi written here,
## which shares no code with the taped version.

## ---- local reference implementations --------------------------------

# scaled forward algorithm; lpmat is T x K of log emission densities
ref_forward <- function(lpmat, Gof, rows, delta) {
  a <- delta * exp(lpmat[rows[1L], ])
  s <- sum(a); llk <- log(s); a <- a / s
  for (k in seq_along(rows)[-1L]) {
    a <- as.vector(a %*% Gof(rows[k - 1L])) * exp(lpmat[rows[k], ])
    s <- sum(a); llk <- llk + log(s); a <- a / s
  }
  llk
}

ref_tpm <- function(eta_row) {
  e <- exp(c(0, eta_row))
  e / sum(e)
}

ref_seq_rows <- function(g, t) {
  lapply(split(seq_along(g), g), function(r) r[order(t[r])])
}

# forward-backward and Viterbi, one sequence at a time, in logs
ref_decode <- function(lpmat, lGof, rows, ldelta) {
  Tl <- length(rows)
  K <- ncol(lpmat)
  lse <- function(z) { m <- max(z); m + log(sum(exp(z - m))) }
  A <- matrix(0, Tl, K); B <- matrix(0, Tl, K)
  Dl <- matrix(0, Tl, K); Ptr <- matrix(0L, Tl, K)
  A[1L, ] <- ldelta + lpmat[rows[1L], ]
  Dl[1L, ] <- A[1L, ]
  for (s in seq_len(Tl)[-1L]) {
    lG <- lGof(rows[s - 1L])
    for (j in seq_len(K)) {
      A[s, j] <- lse(lG[, j] + A[s - 1L, ])
      z <- lG[, j] + Dl[s - 1L, ]
      Ptr[s, j] <- which.max(z)
      Dl[s, j] <- max(z)
    }
    A[s, ] <- A[s, ] + lpmat[rows[s], ]
    Dl[s, ] <- Dl[s, ] + lpmat[rows[s], ]
  }
  if (Tl > 1L) {
    for (s in rev(seq_len(Tl - 1L))) {
      lG <- lGof(rows[s])
      for (i in seq_len(K)) {
        B[s, i] <- lse(lG[i, ] + lpmat[rows[s + 1L], ] + B[s + 1L, ])
      }
    }
  }
  L <- A + B
  P <- exp(L - apply(L, 1, max))
  path <- integer(Tl)
  path[Tl] <- which.max(Dl[Tl, ])
  if (Tl > 1L) {
    for (s in rev(seq_len(Tl - 1L))) path[s] <- Ptr[s + 1L, path[s + 1L]]
  }
  list(P = P / rowSums(P), path = path,
       llk = lse(A[Tl, ]))
}

sim_hmm <- function(N, Tl, G, mu, sigma, seed) {
  set.seed(seed)
  K <- nrow(G)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- sample.int(K, 1L)
    for (t in seq_len(Tl)[-1L]) {
      s[t] <- sample.int(K, 1L, prob = G[s[t - 1L], ])
    }
    data.frame(id = id, t = seq_len(Tl), state = s,
               y = stats::rnorm(Tl, mu[s], sigma[s]))
  }))
}

# the state-suffixed intercepts of a fit, as a flat named vector
icpt <- function(fit) {
  e <- unlist(fixef(fit))
  stats::setNames(unname(e), sub("\\.\\(Intercept\\)$", "", names(e)))
}

G2 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2, 2, byrow = TRUE)

## ---- stage 1: constant transitions -----------------------------------

test_that("a gaussian HMM matches an independent numeric forward", {
  dd <- sim_hmm(25, 24, G2, c(0, 3), c(0.6, 0.6), 4001)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated"), data = dd)
  e <- icpt(fit)
  mu <- c(e[["mu1"]], e[["mu2"]])
  sg <- exp(c(e[["sigma1"]], e[["sigma2"]]))
  G <- rbind(ref_tpm(e[["tr12"]]), ref_tpm(e[["tr22"]]))
  dl <- ref_tpm(fit$estimates[["hmm_ldel"]])
  lpm <- vapply(seq_len(2), function(k)
    stats::dnorm(dd$y, mu[k], sg[k], log = TRUE), numeric(nrow(dd)))
  ll <- sum(vapply(ref_seq_rows(dd$id, dd$t),
                   function(r) ref_forward(lpm, function(z) G, r, dl),
                   numeric(1)))
  expect_equal(ll, as.numeric(logLik(fit)), tolerance = 1e-10)
  expect_equal(attr(logLik(fit), "df"), 7L)
})

test_that("the family= and + spellings give the same fit", {
  dd <- sim_hmm(12, 15, G2, c(0, 3), c(0.6, 0.6), 4002)
  fa <- frm(bf(y ~ 1),
            family = hmm(K = 2, gaussian(), time = t, group = id),
            data = dd)
  fp <- frm(bf(y ~ 1) + hmm(K = 2, gaussian(), time = t, group = id),
            data = dd)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fp)))
  expect_equal(unlist(fixef(fa)), unlist(fixef(fp)))
  # a one-sided formula is accepted for time/group as well
  ff <- frm(bf(y ~ 1),
            family = hmm(K = 2, gaussian(), time = ~t, group = ~id),
            data = dd)
  expect_equal(as.numeric(logLik(ff)), as.numeric(logLik(fa)))
})

test_that("a gaussian HMM agrees with depmixS4", {
  skip_if_not_installed("depmixS4")
  dd <- sim_hmm(25, 24, G2, c(0, 3), c(0.6, 0.6), 4001)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated"), data = dd)
  dm <- depmixS4::depmix(y ~ 1, data = dd, nstates = 2,
                         ntimes = as.integer(table(dd$id)))
  set.seed(11)
  best <- -Inf
  for (i in 1:4) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE, emcontrol = depmixS4::em.control(
        random.start = TRUE, tol = 1e-12, maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      best <- max(best, as.numeric(depmixS4::logLik(ff)))
    }
  }
  expect_equal(best, as.numeric(logLik(fit)), tolerance = 1e-6)
})

test_that("a T = 1 grouping is refused; fixed transitions collapse it to mixture()", {
  set.seed(4003)
  n <- 300
  z <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  dd <- data.frame(y = stats::rnorm(n, c(0, 3)[z], c(0.9, 0.6)[z]),
                   id = seq_len(n), t = 1L)
  expect_error(
    frm(bf(y ~ 1),
        family = hmm(K = 2, gaussian(), time = t, group = id,
                     init = "estimated"), data = dd),
    "every sequence has length 1")
  # held at constants the transition block is identified again, and the
  # model IS a two-component mixture: probe F2 of dev/hmm-feasibility.md
  fh <- frm(bf(y ~ 1, tr12 = 0, tr22 = 0),
            family = hmm(K = 2, gaussian(), time = t, group = id,
                         init = "estimated"), data = dd)
  fm <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
  expect_equal(as.numeric(logLik(fh)), as.numeric(logLik(fm)),
               tolerance = 1e-10)
  eh <- icpt(fh); em <- icpt(fm)
  for (nm in c("mu1", "mu2", "sigma1", "sigma2")) {
    expect_equal(eh[[nm]], em[[nm]], tolerance = 1e-6)
  }
  # the free initial distribution IS the mixing weight
  expect_equal(1 / (1 + exp(fh$estimates[["hmm_ldel"]])),
               unname(stats::plogis(em[["theta1"]])), tolerance = 1e-6)
  expect_equal(attr(logLik(fh), "df"), attr(logLik(fm), "df"))
})

## ---- stage 2: covariate transitions and the init modes ---------------

test_that("covariate transitions match a numeric forward and depmixS4", {
  dd <- sim_hmm(25, 20, G2, c(0, 3), c(0.7, 0.7), 4004)
  set.seed(4005)
  dd$x <- stats::rnorm(nrow(dd))
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated", trans = ~x), data = dd)
  e <- unlist(fixef(fit))
  mu <- c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]])
  sg <- exp(c(e[["sigma1.(Intercept)"]], e[["sigma2.(Intercept)"]]))
  # the covariate at row r drives the step OUT of r (depmixS4's rule)
  Gof <- function(r) {
    rbind(ref_tpm(e[["tr12.(Intercept)"]] + e[["tr12.x"]] * dd$x[r]),
          ref_tpm(e[["tr22.(Intercept)"]] + e[["tr22.x"]] * dd$x[r]))
  }
  lpm <- vapply(seq_len(2), function(k)
    stats::dnorm(dd$y, mu[k], sg[k], log = TRUE), numeric(nrow(dd)))
  dl <- ref_tpm(fit$estimates[["hmm_ldel"]])
  ll <- sum(vapply(ref_seq_rows(dd$id, dd$t),
                   function(r) ref_forward(lpm, Gof, r, dl), numeric(1)))
  expect_equal(ll, as.numeric(logLik(fit)), tolerance = 1e-10)

  skip_if_not_installed("depmixS4")
  dm <- depmixS4::depmix(y ~ 1, data = dd, nstates = 2, transition = ~x,
                         ntimes = as.integer(table(dd$id)))
  set.seed(12)
  best <- -Inf; bp <- NULL
  for (i in 1:5) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE, emcontrol = depmixS4::em.control(
        random.start = TRUE, tol = 1e-12, maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      v <- as.numeric(depmixS4::logLik(ff))
      if (v > best) { best <- v; bp <- depmixS4::getpars(ff) }
    }
  }
  expect_equal(best, as.numeric(logLik(fit)), tolerance = 1e-6)
  # the fitted transition matrix at x = 0, up to the state relabelling
  # the two optimizers may disagree about
  # depmixS4 lays its transition parameters out coefficient by
  # coefficient, each preceded by the zero of the state-1 reference
  # cell: (0, int_row1, 0, slope_row1, 0, int_row2, 0, slope_row2)
  ints <- bp[3:10][c(2, 6)]
  dmG <- rbind(ref_tpm(ints[1L]), ref_tpm(ints[2L]))
  Gx0 <- rbind(ref_tpm(e[["tr12.(Intercept)"]]),
               ref_tpm(e[["tr22.(Intercept)"]]))
  if (sum(abs(dmG - Gx0)) > sum(abs(dmG - Gx0[2:1, 2:1]))) {
    Gx0 <- Gx0[2:1, 2:1]
  }
  expect_equal(unname(dmG), unname(Gx0), tolerance = 1e-3)
})

test_that("a transition covariate at t drives the step from t to t + 1", {
  # the depmixS4 convention (probe B1). Under it the covariate on a
  # sequence's LAST row is never read; under the other reading it would
  # be the FIRST row's that goes unused. Putting the only nonzero value
  # on the last row therefore has to leave the likelihood exactly where
  # a zero covariate does.
  dd <- sim_hmm(15, 12, G2, c(0, 3), c(0.6, 0.6), 4027)
  dd$x_last <- as.numeric(dd$t == 12L) * 3
  dd$x_first <- as.numeric(dd$t == 1L) * 3
  mk <- function(f) {
    suppressWarnings(
      frm(bf(y ~ 1), data = dd,
          family = hmm(K = 2, gaussian(), time = t, group = id,
                       init = "estimated", trans = f)))
  }
  ll0 <- as.numeric(logLik(mk(~1)))
  # the last row's covariate is never read, so its slope is a flat
  # direction and the maximum is exactly the intercept-only one
  expect_equal(as.numeric(logLik(mk(~x_last))), ll0, tolerance = 1e-5)
  # the first row's is read, and buys something
  expect_gt(as.numeric(logLik(mk(~x_first))), ll0 + 1e-3)
})

test_that("the stationary initial distribution solves on the tape", {
  dd <- sim_hmm(20, 25, G2, c(0, 3), c(0.6, 0.6), 4006)
  fst <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "stationary"), data = dd)
  fes <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated"), data = dd)
  fun <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "uniform"), data = dd)
  # one free logit more, and never a worse fit
  expect_equal(attr(logLik(fst), "df") + 1L, attr(logLik(fes), "df"))
  expect_equal(attr(logLik(fun), "df"), attr(logLik(fst), "df"))
  expect_gte(as.numeric(logLik(fes)), as.numeric(logLik(fst)) - 1e-6)

  e <- icpt(fst)
  G <- rbind(ref_tpm(e[["tr12"]]), ref_tpm(e[["tr22"]]))
  dl <- as.vector(solve(t(diag(2) - G) + 1, rep(1, 2)))
  lpm <- vapply(seq_len(2), function(k)
    stats::dnorm(dd$y, c(e[["mu1"]], e[["mu2"]])[k],
                 exp(c(e[["sigma1"]], e[["sigma2"]]))[k], log = TRUE),
    numeric(nrow(dd)))
  ll <- sum(vapply(ref_seq_rows(dd$id, dd$t),
                   function(r) ref_forward(lpm, function(z) G, r, dl),
                   numeric(1)))
  expect_equal(ll, as.numeric(logLik(fst)), tolerance = 1e-10)
})

test_that("init = 'stationary' is refused with varying transitions", {
  dd <- sim_hmm(8, 12, G2, c(0, 3), c(0.6, 0.6), 4007)
  set.seed(4008); dd$x <- stats::rnorm(nrow(dd))
  expect_error(
    frm(bf(y ~ 1),
        family = hmm(K = 2, gaussian(), time = t, group = id,
                     init = "stationary", trans = ~x), data = dd),
    "needs a constant transition matrix")
  # ... and a single overridden cell is enough to trigger it
  expect_error(
    frm(bf(y ~ 1, tr12 ~ x),
        family = hmm(K = 2, gaussian(), time = t, group = id,
                     init = "stationary"), data = dd),
    "needs a constant transition matrix")
})

## ---- stage 3: other emission families ---------------------------------

test_that("poisson emissions agree with depmixS4", {
  set.seed(4009)
  dd <- sim_hmm(25, 22, G2, c(0, 0), c(1, 1), 4009)
  dd$y <- stats::rpois(nrow(dd), c(1.5, 9)[dd$state])
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, poisson(), time = t, group = id,
                          init = "estimated"), data = dd)
  e <- icpt(fit)
  lam <- exp(c(e[["mu1"]], e[["mu2"]]))
  expect_equal(sort(lam), c(1.5, 9), tolerance = 0.4)
  lpm <- vapply(lam, function(l) stats::dpois(dd$y, l, log = TRUE),
                numeric(nrow(dd)))
  G <- rbind(ref_tpm(e[["tr12"]]), ref_tpm(e[["tr22"]]))
  dl <- ref_tpm(fit$estimates[["hmm_ldel"]])
  ll <- sum(vapply(ref_seq_rows(dd$id, dd$t),
                   function(r) ref_forward(lpm, function(z) G, r, dl),
                   numeric(1)))
  expect_equal(ll, as.numeric(logLik(fit)), tolerance = 1e-9)

  skip_if_not_installed("depmixS4")
  dm <- depmixS4::depmix(y ~ 1, data = dd, nstates = 2,
                         family = stats::poisson(),
                         ntimes = as.integer(table(dd$id)))
  set.seed(13)
  best <- -Inf
  for (i in 1:4) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE, emcontrol = depmixS4::em.control(
        random.start = TRUE, tol = 1e-12, maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      best <- max(best, as.numeric(depmixS4::logLik(ff)))
    }
  }
  expect_equal(best, as.numeric(logLik(fit)), tolerance = 1e-6)
})

test_that("categorical (multinomial) emissions agree with depmixS4", {
  skip_if_not_installed("depmixS4")
  dd <- sim_hmm(35, 15, G2, c(0, 0), c(1, 1), 4010)
  P1 <- c(0.55, 0.25, 0.12, 0.08)
  P2 <- c(0.08, 0.15, 0.32, 0.45)
  set.seed(4011)
  cats <- vapply(dd$state, function(s)
    sample.int(4L, 1L, prob = if (s == 1L) P1 else P2), integer(1))
  Y <- matrix(0L, nrow(dd), 4L)
  Y[cbind(seq_len(nrow(dd)), cats)] <- 1L
  dd$Y <- Y
  dd$cf <- factor(cats)
  # a rare category inside a state drives its emission logit to the
  # boundary, which an unpenalized multinomial logit will always do and
  # which the optimizer reports as singular convergence; the optimum
  # itself is depmixS4's (probe E of dev/hmm-feasibility.md)
  fit <- suppressWarnings(
    frm(bf(Y ~ 1),
        family = hmm(K = 2, multinomial(K = 4), time = t, group = id,
                     init = "estimated"), data = dd))
  dm <- depmixS4::depmix(cf ~ 1, data = dd, nstates = 2,
                         family = depmixS4::multinomial("identity"),
                         ntimes = as.integer(table(dd$id)))
  set.seed(14)
  best <- -Inf; bp <- NULL
  for (i in 1:5) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE, emcontrol = depmixS4::em.control(
        random.start = TRUE, tol = 1e-12, maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      v <- as.numeric(depmixS4::logLik(ff))
      if (v > best) { best <- v; bp <- depmixS4::getpars(ff) }
    }
  }
  expect_equal(best, as.numeric(logLik(fit)), tolerance = 1e-6)
  e <- icpt(fit)
  pk <- function(k) {
    eta <- c(0, vapply(2:4, function(j) e[[paste0("mu", j, k)]],
                       numeric(1)))
    exp(eta) / sum(exp(eta))
  }
  # depmixS4's last 8 parameters are the two states' category
  # probabilities, in its own state order
  dmp <- matrix(bp[7:14], 2, 4, byrow = TRUE)
  frmp <- rbind(pk(1), pk(2))
  if (sum(abs(dmp - frmp)) > sum(abs(dmp - frmp[2:1, ]))) {
    frmp <- frmp[2:1, ]
  }
  expect_equal(unname(dmp), unname(frmp), tolerance = 1e-4)
})

test_that("categorical emissions take transition covariates too", {
  dd <- sim_hmm(25, 15, G2, c(0, 0), c(1, 1), 4012)
  set.seed(4013)
  dd$x <- stats::rnorm(nrow(dd))
  cats <- vapply(dd$state, function(s)
    sample.int(3L, 1L,
               prob = if (s == 1L) c(0.6, 0.3, 0.1) else c(0.1, 0.3, 0.6)),
    integer(1))
  Y <- matrix(0L, nrow(dd), 3L)
  Y[cbind(seq_len(nrow(dd)), cats)] <- 1L
  dd$Y <- Y
  fit <- suppressWarnings(
    frm(bf(Y ~ 1),
        family = hmm(K = 2, multinomial(K = 3), time = t, group = id,
                     init = "estimated", trans = ~x), data = dd))
  expect_true(is.finite(as.numeric(logLik(fit))))
  # every transition cell picked up the default predictor
  expect_true(all(c("tr12.x", "tr22.x") %in% names(unlist(fixef(fit)))))
})

## ---- stage 4: random effects -------------------------------------------

test_that("random effects in a state mean compose with the forward pass", {
  set.seed(4014)
  N <- 22L; Tg <- 25L
  b <- stats::rnorm(N, 0, 0.8)
  dd <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tg); s[1L] <- 1L
    for (t in seq_len(Tg)[-1L]) {
      s[t] <- sample.int(2L, 1L, prob = G2[s[t - 1L], ])
    }
    data.frame(id = g, t = seq_len(Tg), state = s,
               y = stats::rnorm(Tg, c(0, 3)[s] + b[g] * (s == 1L), 0.6))
  }))
  dd$gf <- factor(dd$id)
  f0 <- frm(bf(y ~ 1),
            family = hmm(K = 2, gaussian(), time = t, group = id,
                         init = "stationary"), data = dd)
  f1 <- frm(bf(y ~ 1 + (1 | gf)),
            family = hmm(K = 2, gaussian(), time = t, group = id,
                         init = "stationary"), data = dd)
  expect_gt(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  # the MAIN formula reaches every state's mean, so `(1 | gf)` there is
  # one random-effect block per state - two variance components, not one
  expect_equal(attr(logLik(f1), "df"), attr(logLik(f0), "df") + 2L)
  vc <- VarCorr(f1)
  expect_equal(length(vc), 2L)
  expect_true(grepl("^mu1", names(vc)[1L]))
  # one state at a time is spelled the ordinary way
  f2 <- frm(bf(y ~ 1, mu2 ~ 1 + (1 | gf)),
            family = hmm(K = 2, gaussian(), time = t, group = id,
                         init = "stationary"), data = dd)
  expect_equal(length(VarCorr(f2)), 1L)
  expect_true(grepl("^mu2", names(VarCorr(f2))[1L]))
  # decoding still runs, conditional on the random-effect modes
  expect_equal(nrow(hmm_probs(f1)), nrow(dd))
  expect_equal(unname(rowSums(hmm_probs(f1))), rep(1, nrow(dd)),
               tolerance = 1e-12)
})

test_that("hmm agrees with hmmTMB on the stationary fixed-effect model", {
  skip_if_not_installed("hmmTMB")
  dd <- sim_hmm(20, 28, matrix(c(0.85, 0.15, 0.2, 0.8), 2, 2, byrow = TRUE),
                c(0, 3), c(0.6, 0.6), 4015)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "stationary"), data = dd)
  e <- icpt(fit)
  mu <- c(e[["mu1"]], e[["mu2"]])
  sg <- exp(c(e[["sigma1"]], e[["sigma2"]]))
  G <- rbind(ref_tpm(e[["tr12"]]), ref_tpm(e[["tr22"]]))
  # hmmTMB reads a column called `state` as KNOWN states, silently, and
  # then maximizes the complete-data likelihood instead (probe D3)
  dh <- data.frame(ID = dd$id, t = dd$t, y = dd$y)
  hid <- hmmTMB::MarkovChain$new(data = dh, n_states = 2,
                                 initial_state = "stationary")
  # hmmTMB validates rowSums(tpm) == 1 exactly; renormalize away the
  # last-ulp float error that mac BLAS leaves
  hid$update_tpm(G / rowSums(G))
  obs <- hmmTMB::Observation$new(
    data = dh, n_states = 2, dists = list(y = "norm"),
    par = list(y = list(mean = mu, sd = sg)))
  hm <- hmmTMB::HMM$new(obs = obs, hid = hid)
  hm$fit(silent = TRUE)
  expect_equal(hm$llk(), as.numeric(logLik(fit)), tolerance = 1e-8)
  hp <- hm$par()
  expect_equal(as.numeric(hp$obspar["y.mean", , 1]), mu, tolerance = 1e-5)
  expect_equal(as.numeric(hp$obspar["y.sd", , 1]), sg, tolerance = 1e-5)
  expect_equal(as.numeric(hp$tpm[, , 1]), as.numeric(G), tolerance = 1e-5)
})

## ---- stage 5: post-processing ------------------------------------------

test_that("hmm_probs() and hmm_viterbi() match a per-sequence reference", {
  dd <- sim_hmm(20, 20, G2, c(0, 3), c(0.6, 0.6), 4016)
  set.seed(4017); dd$x <- stats::rnorm(nrow(dd))
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated", trans = ~x), data = dd)
  e <- unlist(fixef(fit))
  mu <- c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]])
  sg <- exp(c(e[["sigma1.(Intercept)"]], e[["sigma2.(Intercept)"]]))
  lGof <- function(r) {
    log(rbind(ref_tpm(e[["tr12.(Intercept)"]] + e[["tr12.x"]] * dd$x[r]),
              ref_tpm(e[["tr22.(Intercept)"]] + e[["tr22.x"]] * dd$x[r])))
  }
  lpm <- vapply(seq_len(2), function(k)
    stats::dnorm(dd$y, mu[k], sg[k], log = TRUE), numeric(nrow(dd)))
  ldl <- log(ref_tpm(fit$estimates[["hmm_ldel"]]))
  Pref <- matrix(NA_real_, nrow(dd), 2)
  vref <- integer(nrow(dd))
  llk <- 0
  for (r in ref_seq_rows(dd$id, dd$t)) {
    d <- ref_decode(lpm, lGof, r, ldl)
    Pref[r, ] <- d$P
    vref[r] <- d$path
    llk <- llk + d$llk
  }
  expect_equal(llk, as.numeric(logLik(fit)), tolerance = 1e-9)
  P <- hmm_probs(fit)
  expect_equal(dim(P), c(nrow(dd), 2L))
  expect_equal(colnames(P), c("state1", "state2"))
  expect_equal(unname(P), Pref, tolerance = 1e-12)
  expect_equal(unname(rowSums(P)), rep(1, nrow(dd)), tolerance = 1e-12)
  expect_identical(hmm_viterbi(fit), vref)
  # a well-separated 2-state gaussian chain decodes almost perfectly
  expect_gt(mean(hmm_viterbi(fit) == dd$state), 0.95)
})

test_that("fitted() is the occupancy-weighted mean, not state 1's", {
  dd <- sim_hmm(20, 20, G2, c(0, 3), c(0.6, 0.6), 4018)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated"), data = dd)
  e <- icpt(fit)
  P <- hmm_probs(fit)
  expect_equal(as.numeric(fitted(fit)),
               as.numeric(P %*% c(e[["mu1"]], e[["mu2"]])),
               tolerance = 1e-10)
  # rung 1's failure mode: a constant fitted value at state 1's mean
  expect_gt(stats::sd(fitted(fit)), 1)
  expect_gt(stats::cor(fitted(fit), dd$y), 0.85)
  expect_equal(as.numeric(predict(fit, type = "response")),
               as.numeric(fitted(fit)))
  expect_equal(as.numeric(residuals(fit)),
               as.numeric(dd$y - fitted(fit)), tolerance = 1e-10)
  # pearson divides by the total (between + within state) variance
  rp <- residuals(fit, type = "pearson")
  expect_equal(length(rp), nrow(dd))
  expect_lt(abs(stats::sd(rp) - 1), 0.25)
  # a state's own predictor stays reachable
  expect_equal(unique(round(as.numeric(predict(fit, dpar = "mu2")), 8)),
               round(e[["mu2"]], 8))
})

test_that("simulate() draws a state path and then emissions", {
  dd <- sim_hmm(20, 20, G2, c(0, 3), c(0.6, 0.6), 4019)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated"), data = dd)
  s <- simulate(fit, nsim = 3, seed = 99)
  expect_equal(dim(s), c(nrow(dd), 3L))
  expect_equal(names(s), paste0("sim_", 1:3))
  # the draws reproduce the marginal spread and, crucially, the
  # PERSISTENCE: an i.i.d. mixture draw would have lag-1 correlation 0
  v <- s[[1L]]
  expect_lt(abs(stats::sd(v) - stats::sd(dd$y)), 0.4)
  lag1 <- stats::cor(v[-1], v[-length(v)])
  expect_gt(lag1, 0.4)
  expect_equal(attr(simulate(fit, nsim = 1, seed = 5), "seed"),
               attr(simulate(fit, nsim = 1, seed = 5), "seed"))
})

## ---- stage 6: guards ---------------------------------------------------

test_that("a missing response is masked, not dropped", {
  dd <- sim_hmm(25, 20, G2, c(0, 3), c(0.6, 0.6), 4020)
  dd$y_na <- dd$y
  dd$y_na[dd$t %in% c(7, 8, 14)] <- NA
  fit <- frm(bf(y_na ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id,
                          init = "estimated"), data = dd)
  # every row is kept: the chain never shortens
  expect_equal(stats::nobs(fit), nrow(dd))
  e <- icpt(fit)
  mu <- c(e[["mu1"]], e[["mu2"]])
  sg <- exp(c(e[["sigma1"]], e[["sigma2"]]))
  ms <- is.na(dd$y_na)
  lpm <- vapply(seq_len(2), function(k) {
    v <- stats::dnorm(ifelse(ms, 0, dd$y_na), mu[k], sg[k], log = TRUE)
    ifelse(ms, 0, v)
  }, numeric(nrow(dd)))
  G <- rbind(ref_tpm(e[["tr12"]]), ref_tpm(e[["tr22"]]))
  dl <- ref_tpm(fit$estimates[["hmm_ldel"]])
  ll <- sum(vapply(ref_seq_rows(dd$id, dd$t),
                   function(r) ref_forward(lpm, function(z) G, r, dl),
                   numeric(1)))
  expect_equal(ll, as.numeric(logLik(fit)), tolerance = 1e-9)
  # no residual where there is no observation, but the state is still
  # inferred from the neighbours
  expect_true(all(is.na(residuals(fit)[ms])))
  expect_true(all(is.finite(hmm_probs(fit)[ms, ])))
  # dropping the rows instead gives a DIFFERENT (biased) fit
  f2 <- frm(bf(y_na ~ 1),
            family = hmm(K = 2, gaussian(), time = t, group = id,
                         init = "estimated"), data = dd[!ms, ])
  expect_gt(abs(icpt(f2)[["tr12"]] - e[["tr12"]]), 1e-3)
})

test_that("hmm() refuses the modes and terms it cannot express", {
  dd <- sim_hmm(10, 12, G2, c(0, 3), c(0.6, 0.6), 4021)
  dd$w <- 1; dd$cc <- 0; dd$sdv <- 1
  dd$gf <- factor(dd$id)
  fam <- function() hmm(K = 2, gaussian(), time = t, group = id)
  expect_error(frm(bf(y | weights(w) ~ 1), family = fam(), data = dd),
               "cannot be combined with weights\\(\\)")
  expect_error(frm(bf(y | cens(cc) ~ 1), family = fam(), data = dd),
               "cannot be combined with cens\\(\\)")
  expect_error(frm(bf(y | trunc(lb = -9) ~ 1), family = fam(), data = dd),
               "cannot be combined with trunc\\(\\)")
  expect_error(frm(bf(y | se(sdv) ~ 1), family = fam(), data = dd),
               "cannot be combined with se\\(\\)")
  expect_error(frm(bf(y ~ 1), family = fam(), data = dd, REML = TRUE),
               "REML = TRUE cannot be combined with hmm")
  expect_error(frm(bf(y ~ 1 + (1 | gf)), family = fam(), data = dd,
                   quadrature = TRUE),
               "quadrature = TRUE cannot be combined with hmm")
  expect_error(frm(bf(y ~ 1), family = fam(), data = dd,
                   control = frmtmb_control(profile = TRUE)),
               "profile = TRUE\\) cannot be combined with hmm")
  expect_error(
    frm(mvbf(bf(y ~ 1) + hmm(K = 2, gaussian(), time = t, group = id),
             bf(x ~ 1) + gaussian()),
        data = transform(dd, x = y)),
    "supports univariate models only")
  # duplicated time points inside a sequence
  expect_error(
    frm(bf(y ~ 1),
        family = hmm(K = 2, gaussian(), time = state, group = id),
        data = dd),
    "time points must be unique within a sequence")
})

test_that("hmm() refuses component families it cannot suffix", {
  expect_error(hmm(2, cumulative()), "ordinal families are not supported")
  expect_error(hmm(2, mixture(gaussian(), gaussian())),
               "cannot be a mixture")
  expect_error(hmm(2, hmm(2, gaussian())), "cannot itself be an hmm")
  expect_error(hmm(1, gaussian()), "at least 2")
  expect_error(hmm(12, gaussian()), "at most 9 states")
  expect_error(hmm(2, gaussian(), trans = y ~ x), "one-sided formula")
})

test_that("post-fit methods that cannot be right refuse instead", {
  dd <- sim_hmm(12, 15, G2, c(0, 3), c(0.6, 0.6), 4022)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  expect_error(predict(fit, type = "response", se.fit = TRUE),
               "se.fit is not supported on the response scale for an hmm")
  expect_error(predict(fit, type = "response", newdata = dd),
               "not available for newdata")
  expect_error(residuals(fit, type = "osa"), "not available for an hmm")
  expect_error(residuals(fit, type = "deviance"), "no per-row likelihood")
  expect_error(conditional_effects(fit), "not available for an hmm")
  # a family with an hmm() lpdf must never hand back a per-row density
  expect_error(fit$spec$responses[[1L]]$family$lpdf(1, list(), list()),
               "per-SEQUENCE forward recursion")
})

test_that("a start on the label-symmetry axis warns", {
  dd <- sim_hmm(12, 15, G2, c(0, 3), c(0.6, 0.6), 4023)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  st <- fit$estimates[["beta"]]
  st[] <- 0
  expect_warning(
    frm(bf(y ~ 1),
        family = hmm(K = 2, gaussian(), time = t, group = id),
        data = dd, start = list(beta = st)),
    "fixed point of the label symmetry")
})

test_that("the default start spreads the state means over the response", {
  dd <- sim_hmm(12, 15, G2, c(0, 3), c(0.6, 0.6), 4024)
  fam <- hmm(K = 3, gaussian(), time = t, group = id)
  qs <- vapply(1:3, function(k) fam$init_dpars[[paste0("mu", k)]](dd$y,
                                                                  list()),
               numeric(1))
  expect_equal(qs, unname(stats::quantile(dd$y, (1:3) / 4)))
  # a log-link family whose quantiles leave the link's range still gets a
  # spread start rather than -Inf
  fp <- hmm(K = 2, poisson(), time = t, group = id)
  y0 <- c(rep(0L, 30), 1:10)
  v <- vapply(1:2, function(k) fp$init_dpars[[paste0("mu", k)]](y0, list()),
              numeric(1))
  expect_true(all(is.finite(log(v))))
  expect_gt(abs(diff(log(v))), 0.5)
})

test_that("time and group default to row order and one sequence", {
  dd <- sim_hmm(1, 90, G2, c(0, 3), c(0.6, 0.6), 4025)
  f1 <- frm(bf(y ~ 1),
            family = hmm(K = 2, gaussian(), time = t, group = id),
            data = dd)
  f2 <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian()), data = dd)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
  # rows out of time order are reordered by `time`, not left as they lie
  sh <- dd[sample(nrow(dd)), ]
  f3 <- frm(bf(y ~ 1),
            family = hmm(K = 2, gaussian(), time = t, group = id),
            data = sh)
  expect_equal(as.numeric(logLik(f3)), as.numeric(logLik(f1)),
               tolerance = 1e-8)
  # ONE sequence is the degenerate shape of the time-sliced decoding
  # passes: every step has a single live sequence, and a K-vector must
  # not be allowed to stand in for a 1 x K matrix
  P <- hmm_probs(f1)
  expect_equal(dim(P), c(nrow(dd), 2L))
  expect_equal(unname(rowSums(P)), rep(1, nrow(dd)), tolerance = 1e-12)
  expect_length(hmm_viterbi(f1), nrow(dd))
  expect_true(all(hmm_viterbi(f1) %in% 1:2))
  expect_length(fitted(f1), nrow(dd))
  expect_equal(dim(simulate(f1, nsim = 1, seed = 3)), c(nrow(dd), 1L))
  # against the per-sequence reference
  e <- icpt(f1)
  lpm <- vapply(seq_len(2), function(k)
    stats::dnorm(dd$y, c(e[["mu1"]], e[["mu2"]])[k],
                 exp(c(e[["sigma1"]], e[["sigma2"]]))[k], log = TRUE),
    numeric(nrow(dd)))
  lG <- log(rbind(ref_tpm(e[["tr12"]]), ref_tpm(e[["tr22"]])))
  G <- exp(lG)
  ldl <- log(as.vector(solve(t(diag(2) - G) + 1, rep(1, 2))))
  ref <- ref_decode(lpm, function(z) lG, seq_len(nrow(dd)), ldl)
  expect_equal(unname(P), ref$P, tolerance = 1e-12)
  expect_identical(hmm_viterbi(f1), ref$path)
})

test_that("frm_sample() runs on an hmm fit", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  dd <- sim_hmm(8, 12, G2, c(0, 3), c(0.6, 0.6), 4026)
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  s <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0)))
  expect_s3_class(s, "frmtmb_draws")
})
