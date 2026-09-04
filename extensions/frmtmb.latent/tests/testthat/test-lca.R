# lca(): the classic latent class measurement model (poLCA's).
#
# The reference is poLCA itself on its own shipped data, which is an
# EM implementation of exactly this likelihood: agreement there is the
# strongest statement available, and it is checked on the
# log-likelihood, the item profiles, the class sizes, the posterior
# membership and (for latent class regression) the gating
# coefficients. A hand-rolled optim() reference pins the simulated
# case to optimizer precision, and a one-item fit is checked against
# the saturated single-categorical likelihood it must collapse to.

# ------------------------------------------------------------- helpers

# simulated K = 2, J = 4 binary-item data with well-separated classes
sim_lca_data <- function(seed = 5, n = 400, p1 = 0.35,
                         pr = rbind(c(0.90, 0.85, 0.20, 0.75),
                                    c(0.20, 0.15, 0.85, 0.25))) {
  set.seed(seed)
  cl <- rbinom(n, 1, p1) + 1
  J <- ncol(pr)
  Y <- matrix(0L, n, J)
  for (j in seq_len(J)) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
  dd <- data.frame(x = rnorm(n))
  dd$Y <- Y
  list(dd = dd, cl = cl)
}

# hand-rolled negative log-likelihood for binary items: one gating
# logit, then one logit per class per item
lca_ref_nll <- function(Y, K = 2) {
  n <- nrow(Y)
  J <- ncol(Y)
  function(p) {
    th <- c(p[seq_len(K - 1L)], 0)
    lw <- th - log(sum(exp(th)))
    A <- matrix(p[-seq_len(K - 1L)], K, J)
    L <- matrix(0, n, K)
    for (k in seq_len(K)) {
      s <- rep(0, n)
      for (j in seq_len(J)) {
        t_ <- ifelse(Y[, j] == 2, stats::plogis(A[k, j], log.p = TRUE),
                     stats::plogis(-A[k, j], log.p = TRUE))
        t_[is.na(t_)] <- 0          # a missing item drops its factor
        s <- s + t_
      }
      L[, k] <- s + lw[k]
    }
    m <- apply(L, 1, max)
    -sum(m + log(rowSums(exp(L - m))))
  }
}

# best label permutation of `a`'s rows against `b`, and the distance
align_rows <- function(a, b) {
  K <- nrow(a)
  perms <- if (K == 2L) list(1:2, 2:1) else {
    list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
         c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  }
  d <- vapply(perms, function(p) max(abs(a[p, , drop = FALSE] - b)),
              numeric(1))
  list(perm = perms[[which.min(d)]], dist = min(d))
}

# poLCA agreement needs an optimizer run to the same depth its EM
# reaches at tol = 1e-12
tight <- function() {
  frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                rel.tol = 1e-14, x.tol = 1e-14))
}

# --------------------------------------------------- the fit itself

test_that("lca() fits the measurement model and reports its parts", {
  s <- sim_lca_data()
  fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd)

  # the gating dpars ARE the model's fixed effects; the item profiles
  # are family extra parameters and take no linear predictor
  expect_named(fixef(fit), "theta1")
  # one extra parameter vector per ITEM, K * (C_j - 1) long. Per item
  # rather than one flat vector is what lets every consumer read the
  # item structure back off the fit instead of caching it.
  expect_equal(fit$frame$extra_names, paste0("pi", 1:4))
  for (nm in fit$frame$extra_names) {
    expect_equal(length(fit$estimates[[nm]]), 2L * 1L)
  }

  pf <- lca_profiles(fit)
  expect_s3_class(pf, "frmtmb_lca_profiles")
  expect_length(pf, 4L)
  for (m in pf) {
    expect_equal(dim(m), c(2L, 2L))
    expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-10)
  }
  cs <- attr(pf, "class_sizes")
  expect_length(cs, 2L)
  expect_equal(sum(cs), 1, tolerance = 1e-10)
  expect_output(print(pf), "Estimated class sizes")

  P <- lca_probs(fit)
  expect_equal(dim(P), c(nrow(s$dd), 2L))
  expect_equal(unname(rowSums(P)), rep(1, nrow(s$dd)), tolerance = 1e-10)
  # separated classes classify sharply
  expect_gt(attr(P, "entropy"), 0.7)
  acc <- mean(max.col(P) == s$cl)
  expect_gt(max(acc, 1 - acc), 0.9)

  # lca_probs() is mixture_probs() with the LCA check and the entropy
  expect_equal(unname(unclass(mixture_probs(fit))),
               unname(matrix(P, nrow(P))), tolerance = 1e-12)
})

test_that("lca() matches a hand-rolled optim() reference exactly", {
  s <- sim_lca_data()
  fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd)
  nll <- lca_ref_nll(s$dd$Y)
  st <- c(0, rep(c(1.5, -1.5), 4))
  op <- stats::optim(st, nll, method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 8000))
  op <- stats::optim(op$par, nll, method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 8000))
  expect_lt(abs(op$value + as.numeric(logLik(fit))), 1e-8)
})

test_that("a covariate on the formula gates class membership", {
  set.seed(11)
  n <- 800
  x <- rnorm(n)
  cl <- 1L + rbinom(n, 1, stats::plogis(-0.8 + 1.2 * x))
  pr <- rbind(c(0.90, 0.85, 0.20, 0.75), c(0.15, 0.10, 0.90, 0.20))
  Y <- matrix(0L, n, 4)
  for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
  dd <- data.frame(x = x)
  dd$Y <- Y

  fit <- frm(bf(Y ~ x), family = lca(K = 2), data = dd)
  cf <- fixef(fit)$theta1
  expect_named(cf, c("(Intercept)", "x"))
  # theta1 is class 1 against class 2; the simulated slope puts high x
  # in class 2, so the fitted slope is negative under that labeling
  # and positive under the swap. Its magnitude is what is estimated.
  expect_gt(abs(cf[["x"]]), 0.8)
  expect_lt(abs(abs(cf[["x"]]) - 1.2), 0.35)

  # a covariate-free fit is nested inside it
  fit0 <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(fit0)))
  expect_lt(anova(fit0, fit)$`Pr(>Chisq)`[2], 0.001)
})

# -------------------------------------------------- poLCA agreement

test_that("lca() reproduces poLCA on the carcinoma data (K = 3)", {
  skip_if_not_installed("poLCA")
  suppressMessages(loadNamespace("poLCA"))
  utils::data("carcinoma", package = "poLCA", envir = environment())

  set.seed(11)
  pl <- poLCA::poLCA(cbind(A, B, C, D, E, F, G) ~ 1, carcinoma,
                     nclass = 3, nrep = 20, verbose = FALSE,
                     maxiter = 20000, tol = 1e-12)

  dd <- carcinoma
  dd$Y <- data.matrix(carcinoma)
  # the 3-class carcinoma solution sits on a boundary (item profiles at
  # 0 and 1), so nlminb reports singular convergence at an optimum poLCA
  # reaches too; the log-likelihood comparison below is the check
  fit <- suppressWarnings(
    frm(bf(Y ~ 1), family = lca(K = 3), data = dd, control = tight()))

  # EM and gradient ascent walk to the same optimum
  expect_lt(abs(as.numeric(logLik(fit)) - pl$llik), 1e-6)

  pf <- lca_profiles(fit)
  al <- align_rows(do.call(cbind, pf), do.call(cbind, pl$probs))
  expect_lt(al$dist, 1e-4)      # item profiles to four decimals
  expect_lt(max(abs(attr(pf, "class_sizes")[al$perm] - pl$P)), 1e-4)
  expect_lt(max(abs(lca_probs(fit)[, al$perm] - pl$posterior)), 1e-4)
})

test_that("lca() reproduces poLCA's latent class regression", {
  skip_if_not_installed("poLCA")
  suppressMessages(loadNamespace("poLCA"))
  utils::data("election", package = "poLCA", envir = environment())
  items <- c("MORALG", "CARESG", "KNOWG", "LEADG", "DISHONG", "INTELG")

  set.seed(31)
  pl <- poLCA::poLCA(
    stats::as.formula(paste("cbind(", paste(items, collapse = ", "),
                            ") ~ PARTY")),
    election, nclass = 2, nrep = 10, verbose = FALSE,
    maxiter = 20000, tol = 1e-12)

  dd <- election
  # the items are factors; data.matrix() takes their level order, which
  # is the documented conversion and poLCA's own
  dd$Y <- data.matrix(election[items])
  fit <- suppressMessages(
    frm(bf(Y ~ PARTY), family = lca(K = 2), data = dd, control = tight()))

  # poLCA drops incomplete cases (na.rm = TRUE) and so does na.omit
  expect_equal(stats::nobs(fit), pl$N)
  expect_lt(abs(as.numeric(logLik(fit)) - pl$llik), 1e-6)

  pf <- lca_profiles(fit)
  al <- align_rows(do.call(cbind, pf), do.call(cbind, pl$probs))
  expect_lt(al$dist, 1e-4)

  # gating coefficients: frmtmb's theta1 is class 1 against class K,
  # poLCA's coeff is its class 2 against its class 1, so the two agree
  # up to the label permutation (a sign flip when K = 2)
  cf <- unname(fixef(fit)$theta1)
  sgn <- if (identical(al$perm, 1:2)) -1 else 1
  expect_lt(max(abs(sgn * cf - as.numeric(pl$coeff))), 1e-4)
})

# ---------------------------------------------- structural properties

test_that("a one-item lca collapses to a single categorical", {
  # J = 1 is not identified - a mixture of categoricals over one item
  # IS one categorical - so the value to check is that the fit reaches
  # the saturated single-item log-likelihood and nothing more
  dd <- data.frame(row = 1:300)
  set.seed(9)
  dd$Y <- matrix(sample(1:3, 300, TRUE, prob = c(0.5, 0.3, 0.2)),
                 ncol = 1)
  fit <- suppressWarnings(
    frm(bf(Y ~ 1), family = lca(K = 2), data = dd, control = tight()))
  tb <- table(dd$Y[, 1])
  sat <- sum(tb * log(tb / sum(tb)))
  expect_lt(abs(as.numeric(logLik(fit)) - sat), 1e-6)
  expect_length(lca_profiles(fit), 1L)
})

test_that("items may carry different numbers of categories", {
  set.seed(13)
  n <- 600
  cl <- rbinom(n, 1, 0.45) + 1
  Y <- cbind(
    1L + rbinom(n, 1, c(0.85, 0.15)[cl]),                 # 2 categories
    sample(3L, n, TRUE),                                  # placeholder
    1L + rbinom(n, 1, c(0.20, 0.80)[cl]))                 # 2 categories
  for (i in seq_len(n)) {
    Y[i, 2] <- sample.int(3L, 1L,
                          prob = if (cl[i] == 1) c(0.7, 0.2, 0.1)
                                 else c(0.1, 0.3, 0.6))
  }
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
  pf <- lca_profiles(fit)
  expect_equal(vapply(pf, ncol, integer(1)), c(item1 = 2L, item2 = 3L,
                                               item3 = 2L))
  # ncat is inferred as the largest observed code per item
  expect_equal(lengths(fit$estimates[fit$frame$extra_names]),
               c(pi1 = 2L, pi2 = 4L, pi3 = 2L))
  # declaring the counts gives the same fit
  fit2 <- frm(bf(Y ~ 1), family = lca(K = 2, ncat = c(2, 3, 2)),
              data = dd)
  expect_equal(as.numeric(logLik(fit2)), as.numeric(logLik(fit)),
               tolerance = 1e-8)
  # a declared count ABOVE what is observed adds a free, near-zero cell
  fit3 <- frm(bf(Y ~ 1), family = lca(K = 2, ncat = c(3, 3, 2)),
              data = dd)
  expect_equal(ncol(lca_profiles(fit3)[[1]]), 3L)
  expect_lt(max(lca_profiles(fit3)[[1]][, 3]), 1e-3)
})

test_that("the starting values are deterministic, so labels are stable", {
  s <- sim_lca_data()
  a <- frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd)
  b <- frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd)
  expect_equal(unclass(lca_profiles(a)), unclass(lca_profiles(b)),
               ignore_attr = TRUE)
  # class 1 is the low-score end of the item-code scale: item 3 is the
  # one whose "yes" marks the second simulated class, so class 1 must
  # be the one with the HIGH probability on it under this data
  pf <- lca_profiles(a)
  expect_gt(pf[[3]][1, 2], pf[[3]][2, 2])

  # a perturbed start is how multimodality is checked; on separated
  # data it returns to the same optimum
  set.seed(4)
  p0 <- a$frame$par_template[a$frame$extra_names]
  c2 <- frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd,
            start = lapply(p0, function(v) v + rnorm(length(v))))
  expect_lt(abs(as.numeric(logLik(c2)) - as.numeric(logLik(a))), 1e-5)
})

test_that("one family object fits two differently shaped data sets", {
  # A family object is a VALUE. Caching the resolved item structure in
  # its environment at fit time made a reused object report the second
  # fit's shape for the first fit's profiles, silently and with wrong
  # probabilities. Nothing may be written to the family: the item
  # structure travels in the parameters.
  fam <- lca(K = 2)

  a <- sim_lca_data(seed = 21, n = 400)            # 4 binary items
  set.seed(22)
  n <- 500
  clb <- rbinom(n, 1, 0.5) + 1
  Yb <- matrix(0L, n, 6)                           # 6 items, one ternary
  for (j in 1:5) {
    Yb[, j] <- 1L + rbinom(n, 1, c(0.85, 0.15)[clb])
  }
  for (i in seq_len(n)) {
    Yb[i, 6] <- sample.int(3L, 1L,
                           prob = if (clb[i] == 1) c(0.7, 0.2, 0.1)
                                  else c(0.1, 0.2, 0.7))
  }
  db <- data.frame(row = seq_len(n))
  db$Y <- Yb

  fa <- frm(bf(Y ~ 1), family = fam, data = a$dd)
  fb <- frm(bf(Y ~ 1), family = fam, data = db)

  # each fit reports its OWN shape, in either read order
  expect_equal(vapply(lca_profiles(fb), ncol, integer(1)),
               c(item1 = 2L, item2 = 2L, item3 = 2L, item4 = 2L,
                 item5 = 2L, item6 = 3L))
  pa <- lca_profiles(fa)
  expect_length(pa, 4L)
  expect_equal(vapply(pa, ncol, integer(1)),
               c(item1 = 2L, item2 = 2L, item3 = 2L, item4 = 2L))
  for (m in pa) {
    expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-10)
  }

  # and the numbers are the ones an unshared family object gives
  solo <- frm(bf(Y ~ 1), family = lca(K = 2), data = a$dd)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(solo)),
               tolerance = 1e-10)
  expect_equal(unclass(pa), unclass(lca_profiles(solo)),
               ignore_attr = TRUE)
  expect_equal(unname(lca_probs(fa)), unname(lca_probs(solo)),
               ignore_attr = TRUE)

  # the same for simulate(), which never sees the response
  s1 <- simulate(fa, nsim = 1, seed = 3)[[1]]
  s2 <- simulate(fb, nsim = 1, seed = 3)[[1]]
  expect_equal(ncol(s1), 4L)
  expect_equal(ncol(s2), 6L)
  expect_true(all(s1 %in% 1:2))
  expect_equal(max(s2[, 6]), 3L)
})

test_that("simulate() draws classes then items, and recovers profiles", {
  s <- sim_lca_data(seed = 7, n = 4000)
  fit <- suppressWarnings(frm(bf(Y ~ 1), family = lca(K = 2),
                              data = s$dd))
  sim <- simulate(fit, nsim = 2, seed = 4)
  expect_equal(ncol(sim), 2L)
  y1 <- sim[[1]]
  expect_true(is.matrix(y1))
  expect_equal(dim(y1), dim(s$dd$Y))
  expect_true(all(y1 %in% 1:2))

  # round trip: refitting the draw recovers the profiles it came from
  d2 <- data.frame(row = seq_len(nrow(y1)))
  d2$Y <- y1
  f2 <- suppressWarnings(frm(bf(Y ~ 1), family = lca(K = 2), data = d2))
  al <- align_rows(do.call(cbind, lca_profiles(f2)),
                   do.call(cbind, lca_profiles(fit)))
  expect_lt(al$dist, 0.05)
})

test_that("missing item responses are masked, not imputed", {
  s <- sim_lca_data()
  dm <- s$dd
  dm$Y[1:40, 1] <- NA
  dm$Y[41:60, 3] <- NA

  # the default drops the incomplete subjects, which is poLCA's default
  f_drop <- suppressMessages(
    frm(bf(Y ~ 1), family = lca(K = 2), data = dm))
  expect_equal(stats::nobs(f_drop), nrow(dm) - 60L)

  # na.rm = FALSE keeps them and only the missing item's factor leaves
  f_keep <- frm(bf(Y ~ 1), family = lca(K = 2, na.rm = FALSE), data = dm)
  expect_equal(stats::nobs(f_keep), nrow(dm))
  expect_equal(dim(lca_probs(f_keep)), c(nrow(dm), 2L))

  nll <- lca_ref_nll(dm$Y)
  op <- stats::optim(c(0, rep(c(1.5, -1.5), 4)), nll, method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 8000))
  op <- stats::optim(op$par, nll, method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 8000))
  expect_lt(abs(op$value + as.numeric(logLik(f_keep))), 1e-6)

  # a subject with nothing observed carries no information
  dbad <- dm
  dbad$Y[5, ] <- NA
  expect_error(frm(bf(Y ~ 1), family = lca(K = 2, na.rm = FALSE),
                   data = dbad),
               "no observed item at all")
})

# ------------------------------------------------------------- guards

test_that("lca() refuses what it does not model", {
  s <- sim_lca_data(n = 200)
  expect_error(lca(K = 1), "at least 2")
  expect_error(lca(K = 2.5), "whole number")
  expect_error(lca(K = 2, ncat = 1), "at least 2")

  expect_error(frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd,
                   REML = TRUE), "REML")
  expect_error(frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd,
                   control = frmtmb_control(profile = TRUE)), "profile")
  expect_error(frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd,
                   quadrature = TRUE), "scalar random intercepts")

  dg <- s$dd
  dg$g <- factor(rep(1:20, length.out = nrow(dg)))
  expect_error(frm(bf(Y ~ (1 | g)), family = lca(K = 2), data = dg),
               "growth-mixture")
  expect_error(frm(bf(Y ~ s(x)), family = lca(K = 2), data = dg),
               "growth-mixture")

  dg$w <- 1
  expect_error(frm(bf(Y | weights(w) ~ 1), family = lca(K = 2),
                   data = dg), "addition terms")

  dg$z <- rnorm(nrow(dg))
  expect_error(frm(mvbf(bf(Y ~ 1) + lca(K = 2), bf(z ~ 1) + gaussian()),
                   data = dg), "multivariate")

  dbad <- s$dd
  dbad$Y[1, 1] <- 1.5
  expect_error(frm(bf(Y ~ 1), family = lca(K = 2), data = dbad),
               "whole-number category codes")

  dconst <- s$dd
  dconst$Y[, 2] <- 1L
  expect_error(frm(bf(Y ~ 1), family = lca(K = 2), data = dconst),
               "fewer than two distinct values")

  dover <- s$dd
  expect_error(frm(bf(Y ~ 1), family = lca(K = 2, ncat = c(2, 2, 2, 2)),
                   data = within(dover, Y[1, 1] <- 3L)),
               "above the declared category count")

  fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = s$dd)
  expect_error(fitted(fit), "no fitted mean")
  expect_error(stats::predict(fit, type = "response"), "no fitted mean")
  expect_error(residuals(fit), "no fitted mean")
  expect_error(residuals(fit, type = "osa"), "whole item response")
  expect_error(lca_probs(1), "fitted model")
  expect_error(lca_profiles(1), "fitted model")
  set.seed(3)
  dg2 <- data.frame(y = rnorm(60))
  gm <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dg2)
  expect_error(lca_probs(gm), "lca\\(\\) family")
  expect_error(lca_profiles(gm), "lca\\(\\) family")
})

test_that("what an lca() fit does support keeps working", {
  s <- sim_lca_data(n = 300)
  fit <- frm(bf(Y ~ x), family = lca(K = 2), data = s$dd)

  # predict() defaults to the gating predictor on the link scale
  p <- stats::predict(fit)
  expect_length(p, nrow(s$dd))
  expect_equal(unname(p), unname(stats::predict(fit, dpar = "theta1")))

  ci <- stats::confint(fit)
  expect_true(all(c("theta1_x", "pi1_1") %in% rownames(ci)))
  expect_s3_class(hypothesis(fit, "theta1_x = 0"), "frmtmb_hypothesis")
  expect_s3_class(frm_allfit(fit), "frmtmb_allfit")
  expect_s3_class(stats::update(fit, . ~ 1), "frmtmb_fit")
})

