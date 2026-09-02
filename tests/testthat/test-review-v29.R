# Confirmed defects from the exotic-model wave (animal models,
# phylogenetic regression, meta-analysis, ordinal cs(), multivariate
# |ID| blocks).
#
# The headline one is F1: the whole point of gr(cov = A) is
# heritability-as-ICC, and the block's standard deviation was invisible
# to hypothesis() and variables(), so the idiom could not be written at
# all. The rest are a false-positive check, a prediction that never left
# the link scale, a p-value for a test that was never run, and a
# delta-method interval that degenerated on a boundary correlation.

skip_on_cran()

# A small relationship matrix with real off-diagonal structure: with A
# the identity the genetic and residual variances are NOT separately
# identified, which is exactly what the OLRE test below turns on. The
# related pairs are DISJOINT, so the matrix is block diagonal and stays
# positive definite for any |rho| < 1.
v29_relmat <- function(ng, pairs = 1L, rho = 0.5) {
  A <- diag(ng)
  for (i in seq_len(pairs)) {
    a <- 2L * i - 1L
    b <- 2L * i
    A[a, b] <- A[b, a] <- rho
  }
  dimnames(A) <- list(paste0("a", seq_len(ng)), paste0("a", seq_len(ng)))
  A
}

# One record per individual: the animal-model shape that used to trip
# the observation-level-random-effect warning.
v29_animal_data <- function(seed, ng = 40, per = 1L, sd_u = 0.8,
                            sd_e = 0.5) {
  set.seed(seed)
  A <- v29_relmat(ng, pairs = 8L)
  u <- drop(crossprod(chol(A), stats::rnorm(ng))) * sd_u
  dd <- data.frame(
    id = factor(rep(rownames(A), each = per), levels = rownames(A)),
    x = stats::rnorm(ng * per)
  )
  dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, sd_e)
  list(dd = dd, A = A)
}


## F1 -------------------------------------------------------------------
# hypothesis() and variables() skipped every known-structure block, so
# sd_<group>__<term> did not exist for gr(cov =), gr(prec =) or
# equalto() and the ICC could not be written.

test_that("gr(cov =) blocks expose sd_/cor_ names to hypothesis()", {
  d <- v29_animal_data(1, ng = 40, per = 2L)
  fit <- frm(bf(y ~ x + (1 | gr(id, cov = A))) + gaussian(), data = d$dd,
             data2 = list(A = d$A))

  expect_true("sd_id__Intercept" %in% variables(fit))
  expect_true("sigma" %in% variables(fit))

  # the name is the block's own within-level sd, which is what VarCorr
  # reports for the block
  vc <- as.data.frame(VarCorr(fit))
  h_sd <- hypothesis(fit, "sd_id__Intercept")
  expect_equal(h_sd$estimate, vc$sdcor[1], tolerance = 1e-10)

  # the headline idiom: heritability as an ICC
  h <- hypothesis(fit,
                  "sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2)")
  expect_equal(nrow(h), 1L)
  expect_true(is.finite(h$estimate) && h$estimate > 0 && h$estimate < 1)
  expect_gt(h$se, 0)
  expect_lt(h$lwr, h$estimate)
  expect_gt(h$upr, h$estimate)

  # hand-computed from the same two quantities
  expect_equal(h$estimate,
               vc$sdcor[1]^2 / (vc$sdcor[1]^2 + sigma(fit)^2),
               tolerance = 1e-10)

  # a bootstrap sanity check on the delta-method interval: the two
  # methods agree on the point estimate and overlap substantially
  hb <- hypothesis(fit,
                   "sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2)",
                   method = "boot", nsim = 40, seed = 7)
  expect_equal(hb$estimate, h$estimate, tolerance = 1e-10)
  expect_lt(abs(hb$se - h$se), 0.5 * h$se + 0.05)
  expect_lt(hb$lwr, h$upr)
  expect_gt(hb$upr, h$lwr)
})

test_that("gr(prec =) and equalto() blocks contribute names too", {
  ng <- 12
  Q <- Matrix::bandSparse(ng, k = c(-1, 0, 1),
                          diagonals = list(rep(-0.4, ng - 1),
                                           rep(1.2, ng),
                                           rep(-0.4, ng - 1)))
  dimnames(Q) <- list(paste0("g", seq_len(ng)), paste0("g", seq_len(ng)))
  Sm <- solve(as.matrix(Q))
  set.seed(12)
  b <- drop(crossprod(chol(Sm), stats::rnorm(ng))) * 0.7
  dd <- data.frame(g = factor(rep(rownames(Sm), each = 8),
                              levels = rownames(Sm)),
                   x = stats::rnorm(ng * 8))
  dd$y <- 1 + 0.5 * dd$x + b[as.integer(dd$g)] +
    stats::rnorm(nrow(dd), 0, 0.5)
  fp <- frm(bf(y ~ x + (1 | gr(g, prec = Q))) + gaussian(), data = dd,
            data2 = list(Q = Q))
  expect_true("sd_g__Intercept" %in% variables(fp))
  expect_equal(hypothesis(fp, "sd_g__Intercept")$estimate,
               as.data.frame(VarCorr(fp))$sdcor[1], tolerance = 1e-10)

  # equalto() estimates nothing: its sds are known constants, so the
  # names are there and their delta-method standard errors are zero
  set.seed(13)
  ng2 <- 30
  V <- matrix(c(1.2, 0.5, 0.5, 0.8), 2, 2)
  bb <- t(chol(V)) %*% matrix(stats::rnorm(2 * ng2), 2)
  ed <- data.frame(y = 1 + as.vector(bb) + stats::rnorm(2 * ng2, 0, 0.5),
                   f = factor(rep(c("a", "b"), ng2)),
                   g = factor(rep(seq_len(ng2), each = 2)))
  fe <- frm(bf(y ~ 1 + equalto(f + 0 | g, V)) + gaussian(), data = ed)
  vv <- variables(fe)
  sd_nms <- grep("^sd_g__", vv, value = TRUE)
  expect_length(sd_nms, 2L)
  he <- hypothesis(fe, sd_nms)
  expect_equal(he$se, c(0, 0), tolerance = 1e-10)
  expect_vector_equal(he$estimate, sqrt(diag(V)), tol = 1e-10)
  # the correlation of a fixed covariance is a constant too
  cor_nms <- grep("^cor_g__", vv, value = TRUE)
  expect_length(cor_nms, 1L)
  expect_equal(hypothesis(fe, cor_nms)$estimate,
               V[1, 2] / sqrt(V[1, 1] * V[2, 2]), tolerance = 1e-10)
})

test_that("blocks whose theta is not a set of sds stay excluded", {
  set.seed(14)
  dd <- data.frame(x = stats::runif(120, 0, 5))
  dd$y <- stats::rnorm(120, sin(dd$x), 0.4)
  fs <- frm(bf(y ~ s(x)) + gaussian(), data = dd)
  # a smooth's theta is an inverse smoothing parameter, not an sd
  expect_false(any(grepl("^sd_", variables(fs))))
  # and confint_varcorr() still reports it under its own label
  expect_true(any(grepl("wiggle", confint_varcorr(fs)$term)))
})


## F2 -------------------------------------------------------------------
# vcov(full = TRUE) refused to report theta under REML, which left no
# delta-method route to a variance component for a REML fit.

test_that("vcov(full = TRUE) carries theta under REML", {
  d <- v29_animal_data(2, ng = 30, per = 3L)
  fr <- frm(bf(y ~ x + (1 | gr(id, cov = A))) + gaussian(), data = d$dd,
            data2 = list(A = d$A), REML = TRUE)

  V <- expect_silent(vcov(fr, full = TRUE))
  # the documented invariant: full = TRUE is labeled like confint rows
  expect_identical(rownames(V), rownames(confint(fr)))
  expect_true("theta_1" %in% rownames(V))
  expect_true(all(diag(V) > 0))
  # and it is the covariance confint()'s Wald interval is built from
  expect_equal(sqrt(diag(V)),
               unname((confint(fr)[, "upr"] - confint(fr)[, "lwr"]) /
                        (2 * stats::qnorm(0.975))),
               tolerance = 1e-8, ignore_attr = TRUE)
  # beta is integrated out of the REML outer problem; vcov() still has it
  expect_false("x" %in% rownames(V))
  expect_true("x" %in% rownames(vcov(fr)))

  # so the ICC is available under REML as well
  h <- hypothesis(fr,
                  "sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2)")
  expect_gt(h$se, 0)
})


## F3 -------------------------------------------------------------------
# The observation-level-random-effect check fired on gr(cov = A) fits
# with one row per level, where the fixed relationship matrix breaks the
# confounding it warns about.

test_that("the OLRE check skips known-structure blocks", {
  d <- v29_animal_data(3, ng = 40, per = 1L)

  expect_silent(
    frm(bf(y ~ x + (1 | gr(id, cov = A))) + gaussian(), data = d$dd,
        data2 = list(A = d$A))
  )
  # gr(prec =) is the same argument on the precision side
  Q <- solve(d$A)
  dimnames(Q) <- dimnames(d$A)
  expect_silent(
    frm(bf(y ~ x + (1 | gr(id, prec = Q))) + gaussian(), data = d$dd,
        data2 = list(Q = Q))
  )
  # but a plain (1 | id) with one row per level is still confounded
  expect_warning(
    frm(bf(y ~ x + (1 | id)) + gaussian(), data = d$dd),
    "own random effect"
  )
})


## F4 -------------------------------------------------------------------
# predict(type = "response") returned the linear predictor for the
# ordinal families; an ordinal response has no mean, so what "response"
# has to give is the category distribution.

v29_ordinal_data <- function(seed, n = 250, tau = c(-0.8, 0.6),
                             beta = 0.9) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  eta <- beta * dd$x
  p <- cbind(stats::plogis(tau[1] - eta),
             stats::plogis(tau[2] - eta) - stats::plogis(tau[1] - eta),
             1 - stats::plogis(tau[2] - eta))
  dd$y <- factor(apply(p, 1L, function(pr) sample(3L, 1L, prob = pr)),
                 levels = 1:3, ordered = TRUE)
  dd
}

test_that("predict(type = 'response') gives ordinal category probabilities", {
  dd <- v29_ordinal_data(41)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)

  nd <- data.frame(x = c(-1, 0, 1))
  P <- predict(fit, newdata = nd, type = "response")
  expect_true(is.matrix(P))
  expect_equal(dim(P), c(3L, 3L))
  expect_equal(colnames(P), levels(dd$y))
  expect_equal(unname(rowSums(P)), rep(1, 3), tolerance = 1e-12)

  # hand-computed: differences of plogis(tau_k - eta), with the
  # thresholds read off the fit's own (tau_1, log increment) storage
  raw <- fit$estimates$tau_raw
  tau <- c(raw[1], raw[1] + exp(raw[2]))
  eta <- predict(fit, newdata = nd, type = "link")
  Fm <- cbind(0, stats::plogis(tau[1] - eta),
              stats::plogis(tau[2] - eta), 1)
  ref <- Fm[, -1L, drop = FALSE] - Fm[, -4L, drop = FALSE]
  expect_vector_equal(as.vector(P), as.vector(ref), tol = 1e-12)

  # the probabilities move monotonically with the predictor, as a
  # cumulative model requires
  expect_true(all(diff(P[, 3]) > 0))
  expect_true(all(diff(P[, 1]) < 0))

  # in sample: one row per observation, still a distribution
  Pin <- predict(fit, type = "response")
  expect_equal(nrow(Pin), nobs(fit))
  expect_equal(unname(rowSums(Pin)), rep(1, nobs(fit)), tolerance = 1e-12)

  # the link scale is untouched, and se.fit lives there
  expect_true(is.numeric(predict(fit, newdata = nd)))
  expect_equal(predict(fit, newdata = nd, dpar = "mu", type = "response"),
               predict(fit, newdata = nd, type = "link"))
  expect_error(predict(fit, type = "response", se.fit = TRUE),
               "not supported")
})

test_that("the other three ordinal families predict distributions too", {
  dd <- v29_ordinal_data(42)
  for (fam in list(sratio(), cratio(), acat())) {
    fit <- frm(bf(y ~ x) + fam, data = dd)
    P <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)),
                 type = "response")
    expect_equal(dim(P), c(3L, 3L), info = fam$family)
    expect_equal(unname(rowSums(P)), rep(1, 3), tolerance = 1e-12,
                 info = fam$family)
    expect_true(all(P > 0 & P < 1), info = fam$family)
  }
  # probit is the other supported link
  fp <- frm(bf(y ~ x) + cumulative("probit"), data = dd)
  Pp <- predict(fp, newdata = data.frame(x = 0), type = "response")
  expect_equal(unname(rowSums(Pp)), 1, tolerance = 1e-12)
})

test_that("cs() terms enter the ordinal predictions and are re-evaluated", {
  dd <- v29_ordinal_data(43)
  fit <- frm(bf(y ~ x + cs(x)) + sratio(), data = dd)

  nd <- data.frame(x = c(-1, 0, 1))
  P <- predict(fit, newdata = nd, type = "response")
  expect_equal(unname(rowSums(P)), rep(1, 3), tolerance = 1e-12)
  # non-degenerate: the defect returned a constant (all-zero) prediction
  expect_gt(min(P), 0)
  expect_gt(stats::sd(P[, 1]), 0.05)

  # the cs() coefficients really move the answer: dropping them changes
  # the distribution at a non-zero x
  P0 <- predict(fit, newdata = data.frame(x = 0), type = "response")
  raw <- fit$estimates$tau_raw
  tau <- c(raw[1], raw[1] + exp(raw[2]))
  h <- stats::plogis(tau)          # at x = 0 the cs offsets vanish
  expect_vector_equal(as.vector(P0),
                      c(h[1], (1 - h[1]) * h[2], (1 - h[1]) * (1 - h[2])),
                      tol = 1e-10)

  # in-sample and newdata routes agree on the training rows, which is
  # the check that the newdata re-evaluation of cs(x) is the same column
  Pin <- predict(fit, type = "response")
  Pnd <- predict(fit, newdata = fit$frame$data_frame[1:5, ],
                 type = "response")
  expect_vector_equal(as.vector(Pin[1:5, ]), as.vector(Pnd), tol = 1e-10)
})


## F5 -------------------------------------------------------------------
# anova() ran pchisq(0, df = 0), which is 0, so two models of the same
# dimension printed "< 2.2e-16 ***" for a test that never happened.

test_that("anova() reports NA, not a p-value, at zero df difference", {
  set.seed(51)
  d <- data.frame(x = stats::rnorm(80), z = stats::rnorm(80))
  d$y <- stats::rnorm(80, 1 + 0.5 * d$x, 1)
  m1 <- frm(bf(y ~ x) + gaussian(), data = d)

  a <- anova(m1, m1)
  expect_equal(a$`Chi Df`[2], 0L)
  expect_equal(a$Chisq[2], 0, tolerance = 1e-8)
  expect_true(is.na(a$`Pr(>Chisq)`[2]))
  expect_false(any(grepl("2.2e-16", utils::capture.output(print(a)))))
  expect_false(any(grepl("\\*", utils::capture.output(print(a)))))

  # a non-nested pair of the same dimension is equally untestable
  m2 <- frm(bf(y ~ z) + gaussian(), data = d)
  expect_true(is.na(anova(m1, m2)$`Pr(>Chisq)`[2]))

  # a real df difference still gets its p-value
  m3 <- frm(bf(y ~ x + z) + gaussian(), data = d)
  a3 <- anova(m1, m3)
  expect_equal(a3$`Chi Df`[2], 1L)
  expect_false(is.na(a3$`Pr(>Chisq)`[2]))
})


## F6 -------------------------------------------------------------------
# confint_varcorr() degenerated on a boundary correlation (the Fisher-z
# clamp flattened the jacobian, so lwr == est == upr == the clamp), and
# VarCorr()/ranef() looked their blocks up by name, so two blocks with
# the same term label collapsed onto the first.

test_that("a boundary correlation reports NA bounds, not a clamp", {
  set.seed(61)
  ng <- 12
  u <- stats::rnorm(ng, 0, 0.8)
  per <- 5
  dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                   x = stats::rnorm(ng * per))
  # both responses driven by the SAME group effect, so the estimated
  # correlation of the merged |ID| block runs to 1
  dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.4)
  dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.3)
  fit <- suppressWarnings(
    frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
             bf(y2 ~ x + (1 | q | id)) + gaussian()), data = dd)
  )
  expect_equal(fit$frame$re_blocks[[1]]$dim, 2L)

  expect_warning(confint_varcorr(fit), "boundary")
  ci <- suppressWarnings(confint_varcorr(fit))
  cr <- ci[ci$type == "cor", ]
  expect_equal(nrow(cr), 1L)
  expect_equal(abs(cr$estimate), 1, tolerance = 1e-4)
  # the defect: lwr == upr == 0.9999, a zero-width interval at the clamp
  expect_true(is.na(cr$lwr))
  expect_true(is.na(cr$upr))
  # the sd rows are unaffected and keep real intervals
  sd_rows <- ci[ci$type == "sd", ]
  expect_equal(nrow(sd_rows), 2L)
  expect_true(all(is.finite(c(sd_rows$lwr, sd_rows$upr))))
  expect_true(all(sd_rows$lwr < sd_rows$estimate))
  expect_true(all(sd_rows$upr > sd_rows$estimate))
})

test_that("an ordinary correlation still gets a two-sided interval", {
  set.seed(62)
  n_g <- 30
  g <- factor(rep(seq_len(n_g), each = 8))
  U <- matrix(stats::rnorm(n_g * 2), n_g) %*%
    chol(matrix(c(0.64, 0.2, 0.2, 0.25), 2))
  dd <- data.frame(g = g, x = stats::rnorm(n_g * 8))
  dd$y <- 1 + 0.5 * dd$x + U[g, 1] + U[g, 2] * dd$x +
    stats::rnorm(n_g * 8, 0, 0.6)
  fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
  ci <- expect_silent(confint_varcorr(fit))
  cr <- ci[ci$type == "cor", ]
  expect_lt(cr$lwr, cr$estimate)
  expect_gt(cr$upr, cr$estimate)
})

test_that("blocks sharing a term label are all reported, not the first twice", {
  # an animal model: an additive genetic term and a permanent
  # environment term on the same individual, both deparsing to "1 | id"
  set.seed(63)
  ng <- 26
  A <- v29_relmat(ng, pairs = 13L, rho = 0.7)
  u <- drop(crossprod(chol(A), stats::rnorm(ng))) * 0.8
  pe <- stats::rnorm(ng, 0, 0.35)
  per <- 6
  dd <- data.frame(id = factor(rep(rownames(A), each = per),
                               levels = rownames(A)),
                   x = stats::rnorm(ng * per))
  dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] + pe[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.5)
  fit <- frm(bf(y ~ x + (1 | gr(id, cov = A)) + (1 | id)) + gaussian(),
             data = dd, data2 = list(A = A))

  vc <- VarCorr(fit)
  expect_length(vc, 2L)
  expect_identical(names(vc), c("1 | id", "1 | id"))
  # the two blocks are different objects, not the first one twice
  expect_false(isTRUE(all.equal(vc[[1]], vc[[2]])))

  # the printed table carries BOTH blocks; the defect printed the first
  # one twice, so the two term rows came out identical
  out <- utils::capture.output(print(vc))
  rows <- grep("Intercept", out, value = TRUE)
  expect_length(rows, 2L)
  expect_false(identical(rows[1], rows[2]))
  expect_equal(sum(grepl("1 | id", out, fixed = TRUE)), 2L)

  vdf <- as.data.frame(vc)
  expect_equal(nrow(vdf), 2L)
  expect_equal(length(unique(vdf$sdcor)), 2L)

  # ranef() used to lose a block outright: out[[label]] <- M overwrites
  re <- ranef(fit)
  expect_length(re, 2L)
  expect_false(isTRUE(all.equal(unname(re[[1]]), unname(re[[2]]))))
  expect_equal(nrow(as.data.frame(re)), 2L * ng)

  # confint_varcorr() reports both blocks as well (the permanent
  # environment component is weakly identified against the genetic one
  # on data this small, which is what the wide-interval warning says)
  expect_equal(nrow(suppressWarnings(confint_varcorr(fit))), 2L)
})


## F1 (follow-up) --------------------------------------------------------
# The |ID| guard in parse_linpred() runs BEFORE the gr() rewrite turns
# cls from "us" into "gr_cov"/"gr_prec", so exactly those two slipped
# through; the merged block was then a us density that never read aux_A,
# and the relationship matrix vanished without a symptom. v0.29 refused
# the construct; v0.32 fits it, by building the merged block as one
# gr_cov/gr_prec Kronecker block of the total merged dimension. The
# equivalence with the long-format spelling lives in test-id-kron.R;
# what stays here is that the matrix is genuinely read and that the
# neighboring constructs are unchanged.

v29_mv_gr_data <- function(seed = 6, ng = 20, per = 6) {
  set.seed(seed)
  A <- v29_relmat(ng, pairs = 3L, rho = 0.8)
  u1 <- drop(crossprod(chol(A), stats::rnorm(ng))) * 0.8
  u2 <- drop(crossprod(chol(A), stats::rnorm(ng))) * 0.5
  dd <- data.frame(id = factor(rep(rownames(A), each = per),
                               levels = rownames(A)),
                   x = stats::rnorm(ng * per))
  dd$y1 <- 1 + 0.5 * dd$x + u1[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.7)
  dd$y2 <- -1 + 0.2 * dd$x + u2[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.5)
  list(dd = dd, A = A)
}

test_that("a SHARED |ID| key on a gr() term keeps the matrix", {
  d <- v29_mv_gr_data()
  A <- d$A
  # across responses
  f <- frm(mvbf(bf(y1 ~ x + (1 | q | gr(id, cov = A))) + gaussian(),
                bf(y2 ~ x + (1 | q | gr(id, cov = A))) + gaussian()),
           data = d$dd, data2 = list(A = A))
  expect_length(f$frame$re_blocks, 1L)
  expect_equal(f$frame$re_blocks[[1]]$covstruct, "gr_cov")
  expect_equal(f$frame$re_blocks[[1]]$dim, 2L)
  # the v0.29 symptom: fitting as a plain us block made the likelihood
  # bit-identical to cov = diag(n). It must not be.
  Imat <- diag(nrow(A))
  dimnames(Imat) <- dimnames(A)
  f_I <- frm(mvbf(bf(y1 ~ x + (1 | q | gr(id, cov = Imat))) + gaussian(),
                  bf(y2 ~ x + (1 | q | gr(id, cov = Imat))) + gaussian()),
             data = d$dd, data2 = list(Imat = Imat))
  expect_gt(abs(as.numeric(logLik(f)) - as.numeric(logLik(f_I))), 0.1)

  # across dpars of one response
  # the sigma component is barely identified here, so the merged
  # correlation runs to the boundary; the block structure is the point
  f_dp <- suppressWarnings(
    frm(bf(y1 ~ x + (1 | q | gr(id, cov = A)),
           sigma ~ (1 | q | gr(id, cov = A))) + gaussian(),
        data = d$dd, data2 = list(A = A)))
  expect_equal(f_dp$frame$re_blocks[[1]]$covstruct, "gr_cov")
  # the precision side is the same construct
  Q <- solve(A)
  dimnames(Q) <- dimnames(A)
  f_q <- frm(mvbf(bf(y1 ~ x + (1 | q | gr(id, prec = Q))) + gaussian(),
                  bf(y2 ~ x + (1 | q | gr(id, prec = Q))) + gaussian()),
             data = d$dd, data2 = list(Q = Q))
  expect_equal(f_q$frame$re_blocks[[1]]$covstruct, "gr_prec")
  expect_equal(as.numeric(logLik(f_q)), as.numeric(logLik(f)),
               tolerance = 1e-6)
  # and the compatibility table says so
  expect_equal(frm_compat("|ID|", "gr_cov")$status, "conditional")
  expect_equal(frm_compat("|ID|", "gr_prec")$status, "conditional")
})

test_that("one |ID| label over two grouping specifications is refused", {
  d <- v29_mv_gr_data()
  A <- d$A
  dd <- d$dd
  # a merged block has room for one structure; the key must not name a
  # plain factor in one formula and a relationship matrix in another
  expect_error(
    frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
             bf(y2 ~ x + (1 | q | gr(id, cov = A))) + gaussian()),
        data = dd, data2 = list(A = A)),
    "more than one grouping specification"
  )
  # cov against prec is the same mistake
  Q <- solve(A)
  dimnames(Q) <- dimnames(A)
  expect_error(
    frm(mvbf(bf(y1 ~ x + (1 | q | gr(id, cov = A))) + gaussian(),
             bf(y2 ~ x + (1 | q | gr(id, prec = Q))) + gaussian()),
        data = dd, data2 = list(A = A, Q = Q)),
    "more than one grouping specification"
  )
  # and so is two different matrices
  B <- diag(nrow(A))
  dimnames(B) <- dimnames(A)
  expect_error(
    frm(mvbf(bf(y1 ~ x + (1 | q | gr(id, cov = A))) + gaussian(),
             bf(y2 ~ x + (1 | q | gr(id, cov = B))) + gaussian()),
        data = dd, data2 = list(A = A, B = B)),
    "more than one grouping specification"
  )
})

test_that("an UNSHARED |ID| key keeps its structure and still fits", {
  d <- v29_mv_gr_data()
  A <- d$A
  f_id <- frm(bf(y1 ~ x + (1 | q | gr(id, cov = A))) + gaussian(),
              data = d$dd, data2 = list(A = A))
  expect_equal(f_id$frame$re_blocks[[1]]$covstruct, "gr_cov")

  # a lone key is a no-op: identical to the same term without it
  f_plain <- frm(bf(y1 ~ x + (1 | gr(id, cov = A))) + gaussian(),
                 data = d$dd, data2 = list(A = A))
  expect_equal(as.numeric(logLik(f_id)), as.numeric(logLik(f_plain)),
               tolerance = 1e-10)
  # and the matrix is genuinely used: swapping in the identity moves the
  # likelihood, which is exactly what the merged-block path failed to do
  Imat <- diag(nrow(A))
  dimnames(Imat) <- dimnames(A)
  f_I <- frm(bf(y1 ~ x + (1 | q | gr(id, cov = Imat))) + gaussian(),
             data = d$dd, data2 = list(Imat = Imat))
  expect_gt(abs(as.numeric(logLik(f_id)) - as.numeric(logLik(f_I))), 0.1)
})

test_that("plain |ID| us terms are unaffected by the gr() merge", {
  d <- v29_mv_gr_data()
  f_us <- frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
                   bf(y2 ~ x + (1 | q | id)) + gaussian()), data = d$dd)
  expect_length(f_us$frame$re_blocks, 1L)
  expect_equal(f_us$frame$re_blocks[[1]]$covstruct, "us")
  expect_equal(f_us$frame$re_blocks[[1]]$dim, 2L)
  expect_equal(frm_compat("|ID|", "us")$status, "works")
})

test_that("the long-format spelling fits through the same path", {
  d <- v29_mv_gr_data()
  A <- d$A
  ld <- data.frame(
    id = rep(d$dd$id, 2), x = rep(d$dd$x, 2),
    trait = factor(rep(c("y1", "y2"), each = nrow(d$dd))),
    value = c(d$dd$y1, d$dd$y2)
  )
  f_long <- frm(bf(value ~ 0 + trait + trait:x +
                     (0 + trait | gr(id, cov = A)), sigma ~ 0 + trait) +
                  gaussian(), data = ld, data2 = list(A = A))
  bk <- f_long$frame$re_blocks[[1]]
  # the d > 1 Kronecker path, which is the one that reads aux_A
  expect_equal(bk$covstruct, "gr_cov")
  expect_equal(bk$dim, 2L)
  expect_false(is.null(bk$aux_kron))
  V <- VarCorr(f_long)[[1]]
  expect_equal(dim(V), c(2L, 2L))
  expect_true(all(diag(V) > 0))
})


## F4 (follow-up) --------------------------------------------------------
# Three internal callers of predict() assumed a vector on the response
# scale.

test_that("get_predict keys ordinal categories with a group column", {
  skip_if_not_installed("marginaleffects")
  dd <- v29_ordinal_data(44)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  nd <- data.frame(x = c(-1, 0, 1))

  gp <- marginaleffects::get_predict(fit, newdata = nd)
  # the defect: 9 rows numbered 1..9, silently misaligned with newdata
  expect_equal(nrow(gp), 9L)
  expect_true("group" %in% names(gp))
  expect_equal(sort(unique(gp$rowid)), 1:3)
  expect_equal(sort(unique(gp$group)), levels(dd$y))
  expect_vector_equal(as.numeric(tapply(gp$estimate, gp$rowid, sum)),
                      rep(1, 3), tol = 1e-12)
  # the values are the prediction matrix, flattened the same way
  P <- predict(fit, newdata = nd, type = "response")
  expect_vector_equal(gp$estimate, as.vector(P), tol = 1e-12)

  # the link scale keeps the plain one-row-per-observation shape
  gl <- marginaleffects::get_predict(fit, newdata = nd, type = "link")
  expect_equal(nrow(gl), 3L)
  expect_false("group" %in% names(gl))
  expect_vector_equal(gl$estimate, predict(fit, newdata = nd), tol = 1e-12)

  # a scalar-response fit is untouched
  set.seed(45)
  gd <- data.frame(x = stats::rnorm(60))
  gd$y <- stats::rnorm(60, 1 + 0.5 * gd$x, 1)
  fg <- frm(bf(y ~ x) + gaussian(), data = gd)
  gg <- marginaleffects::get_predict(fg, newdata = gd[1:4, ])
  expect_equal(nrow(gg), 4L)
  expect_false("group" %in% names(gg))
})

test_that("posterior_linpred stays on the mu predictor for an ordinal fit", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- v29_ordinal_data(46, n = 120)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 2))
  nd <- data.frame(x = c(-1, 0, 1))

  pl <- posterior_linpred(ds, newdata = nd, ndraws = 10)
  expect_equal(dim(pl), c(10L, 3L))
  # transform = TRUE is documented as the mu dpar on its natural scale;
  # it used to route through type = "conditional" into the K-column
  # probability matrix
  plt <- posterior_linpred(ds, newdata = nd, ndraws = 10,
                           transform = TRUE)
  expect_equal(dim(plt), c(10L, 3L))
  # the ordinal mu link is the identity, so the two agree exactly
  expect_vector_equal(as.vector(plt), as.vector(pl), tol = 1e-12)

  # posterior_epred is the one that carries the category distribution,
  # as a draws x observations x categories array (brms's convention)
  ep <- posterior_epred(ds, newdata = nd, ndraws = 10)
  expect_equal(dim(ep), c(10L, 3L, 3L))
  expect_null(dimnames(ep)[[1]])
  expect_equal(dimnames(ep)[[3]], levels(dd$y))
  # each draw's slice is a 3 x 3 matrix of distributions
  for (k in seq_len(dim(ep)[1])) {
    expect_vector_equal(rowSums(ep[k, , ]), rep(1, 3), tol = 1e-12)
  }

  # a non-ordinal fit keeps the documented behavior exactly
  set.seed(47)
  gd <- data.frame(x = stats::rnorm(80))
  gd$y <- stats::rpois(80, exp(0.3 + 0.4 * gd$x))
  fp <- frm(bf(y ~ x) + poisson(), data = gd)
  dp <- suppressWarnings(frm_sample(fp, chains = 1, iter = 400,
                                    refresh = 0, seed = 3))
  lp <- posterior_linpred(dp, newdata = nd, ndraws = 5)
  lpt <- posterior_linpred(dp, newdata = nd, ndraws = 5, transform = TRUE)
  expect_vector_equal(as.vector(lpt), exp(as.vector(lp)), tol = 1e-10)
  epp <- posterior_epred(dp, newdata = nd, ndraws = 5)
  expect_equal(dim(epp), c(5L, 3L))
  expect_null(colnames(epp))
})
