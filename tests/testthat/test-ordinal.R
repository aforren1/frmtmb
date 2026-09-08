sim_ord_data <- function(seed = 91, n = 600, K = 4) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- 1.2 * x
  tau <- c(-1, 0.3, 1.4)
  u <- runif(n)
  p <- plogis(outer(rep(1, n), tau) - eta)
  y <- rowSums(u > cbind(p, 1)) + 1L
  data.frame(y = y, x = x)
}

# raw thresholds -> ordered thresholds
tau_of <- function(fit) {
  raw <- fit$estimates$tau_raw
  cumsum(c(raw[1], exp(raw[-1])))
}

test_that("cumulative logit matches MASS::polr", {
  skip_if_not_installed("MASS")
  dd <- sim_ord_data()
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ref <- MASS::polr(factor(y) ~ x, data = dd, method = "logistic",
                    Hess = TRUE)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-5)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-3)
  expect_vector_equal(tau_of(fit), ref$zeta, tol = 1e-3)
})

test_that("cumulative probit matches MASS::polr", {
  skip_if_not_installed("MASS")
  dd <- sim_ord_data(seed = 92)
  fit <- frm(bf(y ~ x) + cumulative(link = "probit"), data = dd)
  ref <- MASS::polr(factor(y) ~ x, data = dd, method = "probit")
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-3)
})

test_that("ordered factor responses and threshold-only models work", {
  dd <- sim_ord_data(seed = 93, n = 300)
  dd$yf <- factor(dd$y, ordered = TRUE)
  f1 <- frm(bf(yf ~ x) + cumulative(), data = dd)
  f2 <- frm(bf(y ~ x) + cumulative(), data = dd)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-8)

  f0 <- frm(bf(y ~ 1) + cumulative(), data = dd)
  # threshold-only: thresholds are the sample cumulative logits
  p <- cumsum(prop.table(table(dd$y)))[-4]
  expect_vector_equal(tau_of(f0), qlogis(p), tol = 1e-4)
})

test_that("ordinal random-intercept model matches a hand-rolled reference", {
  set.seed(94)
  n <- 800
  g <- factor(rep(1:40, 20))
  x <- rnorm(n)
  eta <- 0.8 * x + rnorm(40, 0, 0.9)[g]
  tau <- c(-0.8, 0.8)
  u <- runif(n)
  p <- plogis(outer(rep(1, n), tau) - eta)
  dd <- data.frame(y = rowSums(u > cbind(p, 1)) + 1L, x = x, g = g)
  fit <- frm(bf(y ~ x + (1 | g)) + cumulative(), data = dd)

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  nll_ref <- function(pp) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    nll <- -sum(RTMB::dnorm(pp$u, 0, exp(pp$lsd), log = TRUE))
    eta <- pp$beta * xv + pp$u[gi]
    t2 <- c(pp$t1, pp$t1 + exp(pp$ld))
    Fv <- function(z) 1 / (1 + exp(-z))
    i1 <- as.numeric(yv == 1); i3 <- as.numeric(yv == 3)
    up <- Fv(t2[pmin(yv, 2)] - eta) * (1 - i3) + i3
    lo <- Fv(t2[pmax(yv - 1, 1)] - eta) * (1 - i1)
    nll - sum(log(up - lo))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(beta = 0, t1 = -0.5, ld = 0, lsd = 0,
                              u = numeric(40)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-5)
})

# --- the ordinal link roster, through the registry --------------------

test_that("the ordinal families take brms's link roster and no more", {
  cdf <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  for (l in c(cdf, "softit")) {
    expect_s3_class(cumulative(l), "frmtmb_family")
  }
  for (l in cdf) {
    expect_s3_class(sratio(l), "frmtmb_family")
    expect_s3_class(cratio(l), "frmtmb_family")
  }
  # brms allows softit for cumulative and not for the sequential pair
  expect_error(sratio("softit"), "softit", fixed = TRUE)
  expect_error(cratio("softit"), "softit", fixed = TRUE)
  # acat off the logit is a DIFFERENT density in brms, not a different
  # link, so it is refused rather than silently fitted
  expect_error(acat("probit"), "logit", fixed = TRUE)
  # and the refusal must give the TRUE reason. probit, cloglog,
  # cauchit and softit all map onto (0, 1) and brms accepts every one
  # of them for acat, so the generic ordinal message, which says the
  # link fails to map onto (0, 1), is false here twice over. It was
  # what this constructor used to say.
  for (l in c("probit", "probit_approx", "cloglog", "cauchit", "softit")) {
    m <- tryCatch(acat(l), error = conditionMessage)
    expect_false(grepl("map onto (0, 1); you gave", m, fixed = TRUE),
                 label = paste("acat(", l, ") blames the link"))
    expect_match(m, "brms accepts", fixed = TRUE)
    expect_match(m, "has not written", fixed = TRUE)
  }
  # the families that DO take those links keep the general message,
  # which is true for them
  expect_match(tryCatch(cumulative("log"), error = conditionMessage),
               "map onto (0, 1)", fixed = TRUE)
  for (f in list(cumulative, sratio, cratio, acat)) {
    expect_error(f("log"), "frmtmb-links", fixed = TRUE)
    expect_error(f("identity"), "frmtmb-links", fixed = TRUE)
  }
})

test_that("the ordinal link is reported, and mu's identity link is not", {
  fam <- cumulative("cloglog")
  expect_identical(frmtmb:::family_link_str(fam), "cdf = cloglog")
  # mu really does carry an identity link: the CDF applies to tau - eta
  expect_identical(fam$links$mu$name, "identity")
})

test_that("the ordinal log-space density matches a closed form", {
  # cloglog has an exact difference: F(a) - F(b) = exp(-e^b) - exp(-e^a),
  # so log of it is logspace_sub(-e^b, -e^a) with nothing to cancel.
  # The generalized robust branch is checked against that, and against
  # the plain difference the old logit-only code would have used.
  lk <- frmtmb:::get_link("cloglog")
  q <- lk$logit_eta
  tau <- c(-0.7, 0.4)
  for (eta in seq(-6, 6, by = 0.5)) {
    a <- tau[2] - eta
    b <- tau[1] - eta
    exact <- RTMB::logspace_sub(-exp(b), -exp(a))
    qa <- q(a)
    qb <- q(b)
    rob <- RTMB::logspace_sub(-qb, -qa) + frmtmb:::log_inv_logit(qa) +
      frmtmb:::log_inv_logit(qb)
    expect_equal(rob, exact, tolerance = 1e-10,
                 label = paste("cloglog robust at eta =", eta))
  }
  # and it stays finite where the plain difference has died. The
  # cloglog saturates at tau - eta LARGE, so that is eta large
  # negative: F(9.4) and F(8.3) are both exactly 1 in double precision
  # and their difference is exactly 0, where the truth is exp(-4024).
  far <- -9
  plain <- log(lk$linkinv(tau[2] - far) - lk$linkinv(tau[1] - far))
  qa <- q(tau[2] - far)
  qb <- q(tau[1] - far)
  rob <- RTMB::logspace_sub(-qb, -qa) + frmtmb:::log_inv_logit(qa) +
    frmtmb:::log_inv_logit(qb)
  expect_false(is.finite(plain))
  expect_true(is.finite(rob))
  expect_equal(rob, RTMB::logspace_sub(-exp(tau[1] - far),
                                       -exp(tau[2] - far)),
               tolerance = 1e-10)
})

test_that("a probit cumulative fit agrees with its own category probs", {
  d <- sim_ord_data()
  for (l in c("probit", "cloglog", "cauchit")) {
    fit <- frm(bf(y ~ x), family = cumulative(l), data = d)
    P <- fitted(fit)
    expect_equal(unname(rowSums(P)), rep(1, nrow(d)), tolerance = 1e-10)
    expect_true(all(P > 0))
    # the slope keeps its sign whichever CDF reads the thresholds
    expect_gt(unlist(fixef(fit))[["mu.x"]], 0)
  }
})

test_that("cratio's robust branch reads its CDF where the density does", {
  # it used to read at tau - eta and lean on 1 - F(-x) = F(x), which the
  # cloglog does not satisfy. The fit is compared with the plain path
  # of the same model, reached through cauchit's absence of a robust
  # field only for shape; here the check is that a cloglog cratio fit
  # agrees with a direct evaluation of its own category probabilities.
  d <- sim_ord_data()
  fit <- frm(bf(y ~ x), family = cratio("cloglog"), data = d)
  P <- fitted(fit)
  expect_equal(unname(rowSums(P)), rep(1, nrow(d)), tolerance = 1e-10)
  # the log-likelihood equals the sum of the log of the chosen cell
  ll <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
  expect_equal(ll, as.numeric(logLik(fit)), tolerance = 1e-6)
})

test_that("the taped ordinal density matches the numeric category probs", {
  # two independent paths: the lpdf takes the generalized log-space
  # form off each link's log odds, and ord_cat_probs() builds the whole
  # category distribution from the plain CDF. They have to agree.
  d <- sim_ord_data(n = 300)
  for (fam in c("cumulative", "sratio", "cratio")) {
    links <- c("probit", "cloglog", "cauchit")
    if (fam == "cumulative") links <- c(links, "softit")
    for (l in links) {
      f <- suppressWarnings(frm(bf(y ~ x), family = do.call(fam, list(l)),
                                data = d))
      P <- fitted(f)
      ll <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
      expect_equal(as.numeric(logLik(f)), ll, tolerance = 1e-8,
                   label = paste(fam, l))
    }
  }
})

test_that("cloglog makes cumulative and sratio the same model", {
  # a textbook identity: on the complementary log-log link the
  # cumulative (proportional hazards) model and the stopping-ratio
  # model are the same likelihood. The two reach it through entirely
  # different code, an interior difference of two CDFs against a
  # product of hazards, so agreement is a check on both.
  d <- sim_ord_data(n = 300)
  a <- frm(bf(y ~ x), family = cumulative("cloglog"), data = d)
  b <- frm(bf(y ~ x), family = sratio("cloglog"), data = d)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-7)
  expect_equal(unlist(fixef(a)$mu), unlist(fixef(b)$mu), tolerance = 1e-5)
  # and it is NOT an accident of the data: the logit pair differs
  al <- frm(bf(y ~ x), family = cumulative("logit"), data = d)
  bl <- frm(bf(y ~ x), family = sratio("logit"), data = d)
  expect_gt(abs(as.numeric(logLik(al)) - as.numeric(logLik(bl))), 1e-3)
})
