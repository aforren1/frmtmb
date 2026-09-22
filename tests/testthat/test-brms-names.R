# brms's names and brms's shapes on a maximum-likelihood fit: items 8 to
# 10 of dev/brms-suite-audit.md section 8, the hypothesis() object of
# finding 7 in dev/brms-methods-tests.md, and brms's positional slots on
# fixef(), ranef(), coef() and VarCorr(). dev/brmsnames-findings.md has
# the before-and-after record.

# Tolerances are relative to the reference value, element by element, so
# a small element is not compared at the scale of a large neighbour.
# `bn_exact` is for two routes through the same floating-point
# arithmetic; `bn_fd` for a central finite difference against an exact
# or a Richardson derivative, whose relative error is O(eps^(2/3)), so
# eps^(1/3) leaves a wide margin without letting a wrong derivative
# through. A reference of exactly zero has no relative scale and must be
# matched exactly.
bn_rel <- function(got, ref, tol) {
  expect_equal(length(got), length(ref))
  expect_true(all(abs(got - ref) <= tol * abs(ref)),
              info = paste("worst relative difference",
                           max(abs(got - ref) / abs(ref))))
}
bn_exact <- function(got, ref) bn_rel(got, ref, 64 * .Machine$double.eps)
bn_fd <- function(got, ref) bn_rel(got, ref, .Machine$double.eps^(1 / 3))

bn_data <- function(seed = 31) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(200), g = factor(rep(1:20, 10)))
  u <- cbind(stats::rnorm(20, 0, 0.7), stats::rnorm(20, 0, 0.35))
  dd$y <- stats::rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] +
                         u[dd$g, 2] * dd$x, 1)
  dd
}

test_that("posterior_summary() on a fit refuses and says to sample", {
  dd <- bn_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  # it used to reach posterior_summary.default and die inside
  # as.matrix() with "is.atomic(x) is not TRUE"
  expect_error(posterior_summary(fit), "needs posterior draws")
  expect_error(posterior_summary(fit), "frm_sample", fixed = TRUE)
  fm <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
  expect_error(posterior_summary(fm), "needs draws")
  # the matrix method is untouched, and is brms's on a 3-D array too
  m <- cbind(a = c(1, 2, 3), b = c(2, 4, 9))
  expect_equal(posterior_summary(m)[, "Estimate"], c(a = 2, b = 5))
  A <- array(seq_len(24), c(4, 3, 2),
             dimnames = list(NULL, c("l1", "l2", "l3"), c("c1", "c2")))
  s3 <- posterior_summary(A)
  expect_identical(dim(s3), c(3L, 4L, 2L))
  expect_identical(dimnames(s3),
                   list(c("l1", "l2", "l3"),
                        c("Estimate", "Est.Error", "Q2.5", "Q97.5"),
                        c("c1", "c2")))
  expect_equal(s3["l2", "Estimate", "c2"], mean(A[, "l2", "c2"]))
})

test_that("variables() uses brms's names, b_ on every coefficient", {
  dd <- bn_data()
  fit <- frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = dd)
  # sigma nobody wrote a formula for is brms's `sigma`, on its natural
  # scale, and not also a b_sigma_Intercept
  expect_identical(variables(fit),
                   c("b_Intercept", "b_x", "sd_g__Intercept", "sd_g__x",
                     "cor_g__Intercept__x", "sigma"))
  bn_exact(hypothesis(fit, "sigma", class = NULL)$hypothesis$Estimate,
           sigma(fit))
  # written out as `sigma ~ 1` it is the coefficient b_sigma_Intercept,
  # on the link scale, and there is no `sigma`, as in brms
  f1 <- frm(bf(y ~ x + (1 + x | g), sigma ~ 1) + gaussian(), data = dd)
  expect_true("b_sigma_Intercept" %in% variables(f1))
  expect_false("sigma" %in% variables(f1))
  expect_false("residual__" %in% names(VarCorr(f1)))
  # a distributional parameter's group-level SD carries its name, as in
  # brms. It used to be a second sd_g__Intercept, which the mu block's
  # name hid, so it was not reachable at all
  fs <- frm(bf(y ~ x + (1 | g), sigma ~ (1 | g)) + gaussian(), data = dd)
  v <- variables(fs)
  expect_true(all(c("sd_g__Intercept", "sd_g__sigma_Intercept") %in% v))
  expect_false(isTRUE(all.equal(
    hypothesis(fs, "sd_g__Intercept", class = NULL)$hypothesis$Estimate,
    hypothesis(fs, "sd_g__sigma_Intercept", class = NULL)$hypothesis$Estimate)))
  # the inverse case: the old spelling names nothing
  expect_false(any(c("Intercept", "x", "sigma_Intercept") %in% v))
})

test_that("VarCorr() is brms's structure, keyed by the grouping factor", {
  dd <- bn_data()
  fit <- frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = dd)
  vc <- VarCorr(fit)
  # brms and lme4 both name the entry "g"; frmtmb named it "1 + x | g"
  expect_named(vc, c("g", "residual__"))
  expect_named(vc$g, c("sd", "cor", "cov"))
  expect_identical(dimnames(vc$g$sd),
                   list(c("Intercept", "x"),
                        c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
  expect_identical(dim(vc$g$cor), c(2L, 4L, 2L))
  # frequentist content, stated in ?VarCorr: the estimate, its
  # delta-method standard error, and Wald quantiles
  V <- varcorr_matrices(fit)[[1L]]
  bn_exact(unname(vc$g$sd[, "Estimate"]), unname(sqrt(diag(V))))
  bn_exact(vc$g$cov["x", "Estimate", "Intercept"], V[2, 1])
  expect_identical(vc$g$cor["Intercept", "Estimate", "Intercept"], 1)
  expect_identical(vc$g$cor["Intercept", "Est.Error", "Intercept"], 0)
  bn_exact(unname(vc$g$sd["x", "Q97.5"]),
           unname(vc$g$sd["x", "Estimate"] +
                    stats::qnorm(0.975) * vc$g$sd["x", "Est.Error"]))
  bn_fd(unname(vc$residual__$sd[, "Estimate"]), sigma(fit))
  # probs is honored on a fit, robust and summary = FALSE are refused
  expect_identical(colnames(VarCorr(fit, probs = c(0.1, 0.9))$g$sd),
                   c("Estimate", "Est.Error", "Q10", "Q90"))
  expect_error(VarCorr(fit, robust = TRUE), "median and MAD")
  expect_error(VarCorr(fit, summary = FALSE), "summary = FALSE")
  # no groups and no scalar residual SD: brms's own refusal
  fp <- frm(bf(y ~ x) + poisson(), data = transform(dd, y = rpois(200, 2)))
  expect_error(VarCorr(fp), "does not contain covariance matrices")
})

test_that("brms's summary slot is refused by name, positionally too", {
  dd <- bn_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  # these used to set flatten and condVar and return the default shape,
  # identical() to the call without the argument
  for (call in list(quote(fixef(fit, FALSE)), quote(ranef(fit, FALSE)),
                    quote(coef(fit, FALSE)),
                    quote(VarCorr(fit, NULL, FALSE)))) {
    expect_error(eval(call), "cannot honor summary = FALSE",
                 info = deparse(call))
  }
  # the defaults, spelled out, are accepted and change nothing
  expect_identical(fixef(fit, TRUE), fixef(fit))
  expect_identical(ranef(fit, TRUE, FALSE, c(0.025, 0.975)), ranef(fit))
  # and the frmtmb arguments that moved past `...` still work by name
  expect_named(fixef(fit, flatten = TRUE), c("(Intercept)", "x",
                                             "sigma_(Intercept)"))
  expect_false(is.null(attr(ranef(fit, condVar = TRUE)$g, "condSD")))
})

test_that("the fit methods' leading slots are brms's", {
  skip_if_not_installed("brms")
  # position by position, as far as brms's own arguments go before its
  # `...`; frmtmb's extras must come after
  first_div <- function(ours, theirs) {
    k <- match("...", theirs) - 1L
    d <- which(ours[seq_len(k)] != theirs[seq_len(k)])
    if (length(d)) d[[1L]] else 0L
  }
  # the inverse case, so this cannot pass by comparing nothing: the
  # signature fixef() shipped with diverges at position 2
  expect_identical(first_div(c("object", "flatten", "..."),
                             names(formals(brms:::fixef.brmsfit))), 2L)
  pairs <- list(
    fixef = c("fixef.frmtmb_fit", "fixef.brmsfit"),
    ranef = c("ranef.frmtmb_fit", "ranef.brmsfit"),
    coef = c("coef.frmtmb_fit", "coef.brmsfit"),
    VarCorr = c("VarCorr.frmtmb_fit", "VarCorr.brmsfit"),
    hypothesis = c("hypothesis.frmtmb_fit", "hypothesis.brmsfit"))
  for (nm in names(pairs)) {
    ours <- names(formals(getFromNamespace(pairs[[nm]][1L], "frmtmb")))
    theirs <- names(formals(getFromNamespace(pairs[[nm]][2L], "brms")))
    expect_identical(first_div(ours, theirs), 0L, info = nm)
  }
})

test_that("hypothesis() returns brms's object, with frequentist content", {
  dd <- bn_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  h <- hypothesis(fit, c("x = 0.5", "x > 0", two = "Intercept < 2"))
  expect_s3_class(h, "frmtmb_hypothesis")
  expect_false(inherits(h, "brmshypothesis"))
  expect_named(h, c("hypothesis", "samples", "prior_samples", "class",
                    "alpha"))
  hs <- h$hypothesis
  expect_named(hs, c("Hypothesis", "Estimate", "Est.Error", "CI.Lower",
                     "CI.Upper", "Evid.Ratio", "Post.Prob", "Star"))
  # brms's labels, and a name on the vector replaces one
  expect_identical(hs$Hypothesis, c("(x)-(0.5) = 0", "(x) > 0", "two"))
  expect_identical(h$class, "b")
  est <- unname(fixef_by_dpar(fit)$mu[["x"]])
  se <- unname(sqrt(vcov(fit)["x", "x"]))
  bn_exact(hs$Estimate[1:2], c(est - 0.5, est))
  bn_fd(hs$Est.Error[2], se)
  # brms's interval: 1 - alpha for "=", 1 - 2 alpha for a directional row
  bn_fd(hs$CI.Upper[1], est - 0.5 + stats::qnorm(0.975) * se)
  bn_fd(hs$CI.Upper[2], est + stats::qnorm(0.95) * se)
  # no posterior, so no evidence ratio
  expect_true(all(is.na(hs$Evid.Ratio)) && all(is.na(hs$Post.Prob)))
  # Star: a directional row stars when its one-sided test rejects
  p <- attr(h, "test")$p
  expect_identical(hs$Star[2], if (p[2] < 0.05) "*" else "")
  bn_fd(p[2], stats::pnorm(-est / se))
  expect_identical(nrow(h$samples), 0L)
  # brms's arguments a fit cannot honor are refused by name
  expect_error(hypothesis(fit, "x > 0", robust = TRUE), "median and MAD")
  expect_error(hypothesis(fit, "Intercept > 0", scope = "ranef",
                          group = "g"), "scope")
  # brms's third slot is class
  expect_identical(hypothesis(fit, "Intercept > 0", "sd", "g")$class,
                   "sd_g")
})

test_that("VarCorr() Est.Error is the delta method, checked independently", {
  skip_if_not_installed("numDeriv")
  # Independent of hypothesis() and of VarCorr()'s own finite difference:
  # numDeriv's Richardson Jacobian of varcorr_matrices(), the covariance
  # from vcov(full = TRUE) selected by row name, and residual__ from the
  # log link by hand (the construction of dev/brmsnames-rev-varcorr-se.R).
  # Seed 33 converges cleanly; at 31 the slope sd is at its boundary and
  # the joint covariance has a negative variance
  dd <- bn_data(33)
  expect_identical(frm(bf(y ~ x + (1 + x | g)) + gaussian(),
                       data = dd)$opt$convergence, 0L)
  fit <- frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = dd)
  vc <- VarCorr(fit)
  V <- vcov(fit, full = TRUE)
  th <- fit$estimates[["theta"]]
  thn <- grep("^theta", rownames(V), value = TRUE)
  f <- function(t) {
    M <- varcorr_matrices(fit, t)[[1L]]
    c(sqrt(diag(M)), stats::cov2cor(M)[2, 1])
  }
  J <- numDeriv::jacobian(f, th)
  se <- sqrt(rowSums((J %*% V[thn, thn, drop = FALSE]) * J))
  bn_fd(unname(vc$g$sd[, "Est.Error"]), se[1:2])
  bn_fd(vc$g$cor["x", "Est.Error", "Intercept"], se[3])
  sn <- grep("sigma", rownames(V), value = TRUE)
  expect_length(sn, 1L)
  bn_fd(unname(vc$residual__$sd[, "Est.Error"]),
        sigma(fit) * sqrt(V[sn, sn]))
})

test_that("hypothesis() reads x:fe as brms does, not as R's `:`", {
  # b_x and b_x:fe differ by less than 1 here, so R's `:` returned b_x,
  # a number, with no error (0.58.0 gives 0.2415 where b_x:fe is 1.0794)
  set.seed(8)
  n <- 400
  d <- data.frame(x = stats::rnorm(n),
                  f = factor(sample(c("a", "e"), n, TRUE)))
  d$y <- 1 + 0.3 * d$x + 0.8 * (d$f == "e") + 1.0 * d$x * (d$f == "e") +
    stats::rnorm(n)
  fit <- frm(bf(y ~ x * f) + gaussian(), data = d)
  fe <- fixef_by_dpar(fit)$mu
  expect_lt(abs(fe[["x"]] - fe[["x:fe"]]), 1)
  expect_true("b_x:fe" %in% variables(fit))
  h <- hypothesis(fit, c("x:fe > 0", "x:fe - x = 0"))$hypothesis
  bn_exact(h$Estimate, c(fe[["x:fe"]], fe[["x:fe"]] - fe[["x"]]))
  expect_identical(h$Hypothesis, c("(x:fe) > 0", "(x:fe-x) = 0"))
  bn_exact(hypothesis(fit, "b_x:fe = 0", class = NULL)$hypothesis$Estimate,
           fe[["x:fe"]])
  # brms's check: under class = "b" a full name is b_b_x:fe
  expect_error(hypothesis(fit, "b_x:fe > 0"),
               "cannot be found in the model: \n'b_b_x:fe'", fixed = TRUE)
  expect_error(hypothesis(fit, "nope > 0"), "b_nope")
  fm <- frm_multiple(bf(y ~ x * f) + gaussian(), data = list(d, d))
  bn_exact(hypothesis(fm, "x:fe > 0")$hypothesis$Estimate, fe[["x:fe"]])
})

test_that("names pass through brms's renaming", {
  set.seed(41)
  n <- 240
  dd <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                   f = factor(rep(c("a b", "c-d", "e"), length.out = n)),
                   gs = factor(paste("lev", rep(1:8, length.out = n))),
                   g = factor(rep(1:6, length.out = n)),
                   h = factor(rep(c("p", "q"), each = 3, length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + stats::rnorm(8)[dd$gs])
  fit <- frm(bf(y ~ x + I(x^2) + f + (1 | gs) + (1 | g:h)) + gaussian(),
             data = dd)
  v <- variables(fit)
  expect_true(all(c("b_IxE2", "b_fcMd", "b_fe", "sd_gs__Intercept",
                    "sd_g:h__Intercept") %in% v))
  expect_named(VarCorr(fit), c("gs", "g:h", "residual__"))
  lab <- brms_par_labels(fit)
  # rename_re_levels(): whitespace in a level is a dot; combine_groups():
  # an interaction level joins its parts with _
  expect_true("r_gs[lev.1,Intercept]" %in% lab)
  expect_true("r_g:h[1_p,Intercept]" %in% lab)
  expect_false(anyDuplicated(lab) > 0L)
  bn_exact(hypothesis(fit, "IxE2 + fcMd = 0")$hypothesis$Estimate,
           sum(fixef_by_dpar(fit)$mu[c("I(x^2)", "fc-d")]))

  # make_stan_names(): a response loses its _ and .
  dd$y_a <- dd$y
  dd$y.b <- stats::rnorm(n)
  fm <- frm(mvbf(bf(y_a ~ x), bf(y.b ~ 1), rescor = TRUE) + gaussian(),
            data = dd)
  expect_identical(variables(fm),
                   c("b_ya_Intercept", "b_ya_x", "b_yb_Intercept",
                     "sigma_ya", "sigma_yb", "rescor__ya__yb"))
  vc <- VarCorr(fm)
  expect_named(vc, "residual__")
  expect_identical(rownames(vc$residual__$sd), c("ya", "yb"))
  expect_identical(dimnames(vc$residual__$cor)[[1L]], c("ya", "yb"))
  bn_exact(vc$residual__$cor["yb", "Estimate", "ya"],
           rescor_matrix(fm)[2, 1])
  # brms has no residual__ once any response predicts sigma
  fp <- frm(mvbf(bf(y_a ~ x, sigma ~ x), bf(y.b ~ 1), rescor = FALSE) +
              gaussian(), data = dd)
  expect_error(VarCorr(fp), "does not contain covariance matrices")
})

test_that("name collisions: refused or suffixed, as brms does each", {
  set.seed(12)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  sigma_z = stats::rnorm(n), Intercept = stats::rnorm(n),
                  g = factor(rep(1:10, 30)))
  d$id <- d$g
  d$y <- 1 + 2 * d$sigma_z + stats::rnorm(n, 0, exp(0.4 * d$z))
  # brms: "Internal renaming led to duplicated names"
  expect_error(frm(bf(y ~ Intercept + x) + gaussian(), data = d),
               "Internal renaming led to duplicated names")
  # brms: the later one takes repair_stanfit()'s __1
  fit <- frm(bf(y ~ sigma_z + (1 | g), sigma ~ z) + gaussian(), data = d)
  v <- variables(fit)
  expect_true(all(c("b_sigma_z", "b_sigma_z__1") %in% v))
  expect_false(anyDuplicated(v) > 0L)
  expect_false(anyDuplicated(brms_par_labels(fit)) > 0L)
  bn_exact(hypothesis(fit, "sigma_z = 0")$hypothesis$Estimate,
           fixef_by_dpar(fit)$mu[["sigma_z"]])
  bn_exact(hypothesis(fit, "sigma_z__1 = 0")$hypothesis$Estimate,
           fixef_by_dpar(fit)$sigma[["z"]])
  # brms: "Duplicated group-level effects are not allowed" when two terms
  # share a coefficient, and the copy of the grouping column is the way
  # to fit both. Not an exact twin such as (1 | g) + (1 | g): brms 2.23.0
  # collapses that one to a single term
  expect_error(frm(bf(y ~ x + (1 | g) + (x | g)) + gaussian(), data = d),
               "Duplicated group-level effects are not allowed")
  f2 <- frm(bf(y ~ 1 + (1 | g) + (1 | id)) + gaussian(), data = d)
  expect_true(all(c("sd_g__Intercept", "sd_id__Intercept") %in%
                    variables(f2)))
})

test_that("frm_simulate(newparams =) takes brms's names only", {
  dd <- bn_data()
  f <- bf(y ~ x + (1 | g)) + gaussian()
  s1 <- frm_simulate(f, dd, newparams = list(b_Intercept = 1, b_x = 0.5,
                                             sigma = 0.7,
                                             sd_g__Intercept = 0.5),
                     nsim = 1, seed = 1)
  expect_length(s1[[1L]], nrow(dd))
  expect_error(frm_simulate(f, dd, newparams = list(Intercept = 1, x = 0.5,
                                                    sigma = 0.7,
                                                    sd_g__Intercept = 0.5)),
               "Intercept -> b_Intercept, x -> b_x", fixed = TRUE)
})

bn_mixture_data <- function() {
  # dev/brmsnames-rev2-mixture.R, data seed 64: 87.75% of rows from
  # component 1, so the share (about 0.88) and its log ratio (about 2)
  # cannot be mistaken for each other. At a share near 0.66 the log
  # ratio is near 0.68 and the defect hides
  set.seed(64)
  n <- 400
  z <- stats::runif(n) < 0.9
  d <- data.frame(x = stats::rnorm(n))
  d$y <- ifelse(z, -2 + 0.3 * d$x + stats::rnorm(n, 0, 0.8),
                2.5 + 0.3 * d$x + stats::rnorm(n, 0, 1.2))
  d
}

test_that("a mixture's theta1, theta2 are brms's simplex, not log ratios", {
  d <- bn_mixture_data()
  fit <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
  v <- variables(fit)
  expect_true(all(c("theta1", "theta2") %in% v))
  expect_false(any(c("b_theta1_Intercept", "theta1_Intercept") %in% v))
  # the log ratio of component 1 against component 2, read by position
  # from the estimate vector, not through any name
  eta <- unname(fit$estimates$betad[[3L]])
  expect_identical(names(fit$estimates$betad)[3L], "theta1_(Intercept)")
  expect_gt(eta, 1.5)
  p1 <- stats::plogis(eta)
  # theta1 alone first: a log ratio read under brms's name gives 1.466
  # here, not 0.377
  h <- hypothesis(fit, "theta1 = 0.5", class = NULL)$hypothesis
  bn_exact(h$Estimate, p1 - 0.5)
  bn_exact(hypothesis(fit, "theta1 + theta2 = 0",
                      class = NULL)$hypothesis$Estimate, 1)
  # the delta method through the share: p (1 - p) times the log ratio's
  # standard error
  se_eta <- sqrt(vcov(fit, full = TRUE)["theta1_(Intercept)",
                                        "theta1_(Intercept)"])
  bn_fd(h$Est.Error[1L], p1 * (1 - p1) * se_eta)

  # three components: theta1..theta3 sum to one, and each is the softmax
  # of the two log ratios against the last component
  set.seed(65)
  d3 <- data.frame(y = c(stats::rnorm(240, -4), stats::rnorm(90, 0),
                         stats::rnorm(70, 4)))
  f3 <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian(), gaussian()),
            data = d3)
  e3 <- unname(f3$estimates$betad[grep("^theta", names(f3$estimates$betad))])
  w <- exp(c(e3, 0)) / sum(exp(c(e3, 0)))
  h3 <- hypothesis(f3, c("theta1 = 0", "theta2 = 0", "theta3 = 0"),
                   class = NULL)$hypothesis
  bn_exact(h3$Estimate, w)
})

test_that("frm_simulate(newparams =) takes brms's mixture simplex", {
  d <- bn_mixture_data()
  f <- bf(y ~ x) + mixture(gaussian(), gaussian())
  np <- list(b_mu1_Intercept = -5, b_mu1_x = 0, b_mu2_Intercept = 5,
             b_mu2_x = 0, sigma1 = 0.5, sigma2 = 0.5,
             theta1 = 0.9, theta2 = 0.1)
  s <- frm_simulate(f, d, newparams = np, nsim = 1, seed = 2)
  share <- mean(s[[1L]] < 0)
  # read as a log ratio, 0.9 would give plogis(0.9) = 0.71
  expect_lt(abs(share - 0.9), abs(share - stats::plogis(0.9)))
  np1 <- np
  np1$theta2 <- NULL
  expect_error(frm_simulate(f, d, newparams = np1),
               "give theta1, theta2, which sum to one", fixed = TRUE)
  np2 <- np
  np2$theta2 <- 0.2
  expect_error(frm_simulate(f, d, newparams = np2), "must sum to one")
})

test_that("brms's names: r_ repeats are suffixed, merged levels refused", {
  set.seed(52)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(paste("lvl", rep(1:12, length.out = n))),
                  f2 = factor(sample(c("u", "v"), n, TRUE)),
                  xo = factor(sample(1:4, n, TRUE), ordered = TRUE))
  d$y <- 1 + 0.4 * d$x + stats::rnorm(12)[as.integer(d$g)] + stats::rnorm(n)
  # dev/brmsnames-rev2-collide.R C4: `lvl 1` and `lvl.1` are both
  # r_gd[lvl.1,Intercept] in brms, which suffixes the later level
  gd <- as.character(d$g)
  i1 <- which(gd == "lvl 1")
  gd[i1[c(TRUE, FALSE)]] <- "lvl.1"
  d$gd <- factor(gd)
  fit <- frm(bf(y ~ x + (1 | gd)) + gaussian(), data = d)
  lab <- brms_par_labels(fit)
  expect_false(anyDuplicated(lab) > 0L)
  r <- grep("^r_gd", lab, value = TRUE)
  expect_identical(r[length(r)], "r_gd[lvl.1,Intercept]__1")
  expect_identical(r[1L], "r_gd[lvl.1,Intercept]")
  expect_identical(levels(d$gd)[c(1L, 13L)], c("lvl 1", "lvl.1"))

  # C1: brms joins the parts of 1_2:3 and 1:2_3 into the one level
  # 1_2_3 and pools them; refused here
  d$gi <- factor(sample(c("1_2", "1"), n, TRUE))
  d$hi <- factor(ifelse(d$gi == "1", "2_3", "3"))
  expect_error(frm(bf(y ~ x + (1 | gi:hi)) + gaussian(), data = d),
               "'1:2_3' and '1_2:3' of group 'gi:hi' are one level in brms",
               fixed = TRUE)

  # C5: brms "Cannot use the same response variable twice"
  d$y_a <- d$y + stats::rnorm(n)
  d$ya <- d$y + stats::rnorm(n)
  expect_error(frm(mvbf(bf(y_a ~ x), bf(ya ~ x), rescor = FALSE) +
                     gaussian(), data = d),
               "Cannot use the same response variable twice")

  # a by-smooth's sds_ drops the ':' that its bs_ keeps, as brms names
  # them; a monotonic scale is bsp_
  fs <- frm(bf(y ~ s(z, by = f2) + mo(xo)) + gaussian(), data = d)
  vs <- variables(fs)
  expect_true(all(c("bs_sz:f2u_1", "sds_szf2u_1", "sds_szf2v_1",
                    "bsp_moxo") %in% vs))
  expect_false(any(c("sds_sz:f2u_1", "b_moxo") %in% vs))
})
