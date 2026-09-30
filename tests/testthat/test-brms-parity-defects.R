# The defect rows of the ported brms suite that lane wt-defects fixed
# (dev/defects-findings.md). Each test pins brms 2.23.0's behavior on
# a model built here, so it runs in the ordinary suite; the ported row
# it answers is named in each test.

tol_eps <- function(ref) 1e3 * .Machine$double.eps * max(1, abs(ref))

test_that("a model-building helper inside a formula is refused by name", {
  # brm:102. brms: "Function 'set_rescor' should not be part of the
  # right-hand side of a formula"; R's "invalid type (list)" before
  set.seed(1)
  dd <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_error(frm(y ~ x + set_rescor(TRUE), data = dd),
               "Function 'set_rescor' should not be part of the right-hand",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(frm(bf(y ~ x + lf(sigma ~ x)), data = dd),
               "Function 'lf' should not be part", fixed = TRUE)
  # the absent case: a user's own function of the same name is a term
  lf <- function(v) v^2
  f <- frm(y ~ lf(x), data = dd)
  expect_identical(names(fixef_by_dpar(f)$mu), c("(Intercept)", "lf(x)"))
})

test_that("categorical() reads numeric codes as categories", {
  # standata:83, :85, :142 and priors:101: "fewer than two categories"
  # of a response with ten
  set.seed(2)
  n <- 120
  dd <- data.frame(x = rnorm(n))
  dd$f <- factor(sample(c("11", "12", "13"), n, replace = TRUE))
  dd$y <- as.numeric(as.character(dd$f))
  ff <- frm(bf(f ~ x), family = categorical(), data = dd)
  fy <- frm(bf(y ~ x), family = categorical(), data = dd)
  expect_named(fixef_by_dpar(fy), c("mu12", "mu13"))
  ll <- as.numeric(logLik(ff))
  expect_lt(abs(as.numeric(logLik(fy)) - ll), tol_eps(ll))
  # brms's bernoulli message for a two-category response
  expect_message(frm(y ~ 1, data = data.frame(y = rep(0:1, 5)),
                     family = categorical(), dry_run = "frame"),
                 "family 'bernoulli' might be a more efficient choice")
  # the prior table names the categories' dpars, as brms's does
  gp <- default_prior(y ~ x, data = dd, family = categorical())
  expect_true(all(c("mu12", "mu13") %in% gp$dpar))
})

test_that("numeric codes of a fitted factor in newdata are its levels", {
  # brmsfit-methods:301, :305, :737; data-helpers:7
  set.seed(3)
  n <- 80
  dd <- data.frame(x = rnorm(n), g = factor(sample(0:1, n, TRUE)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + 0.8 * (dd$g == "1"))
  fit <- frm(y ~ x + g, data = dd)
  nd_f <- data.frame(x = c(0, 1), g = factor(c("0", "1")))
  nd_n <- data.frame(x = c(0, 1), g = c(0, 1))
  expect_identical(fitted(fit, newdata = nd_n), fitted(fit, newdata = nd_f))
  expect_identical(frm_linpred(fit, newdata = nd_n),
                   frm_linpred(fit, newdata = nd_f))
  # a code that is not a level is still a new level
  expect_error(fitted(fit, newdata = data.frame(x = 0, g = 2)),
               "that the fit did not see: '2'", fixed = TRUE)
  # a slice of data whose factor carries a "contrasts" attribute
  # predicts without model.frame()'s "contrasts dropped" warning
  d2 <- dd
  contrasts(d2$g) <- contr.treatment(2)
  fit2 <- frm(y ~ x + g, data = d2)
  expect_no_warning(p <- fitted(fit2, newdata = d2[1:5, ]))
  expect_identical(unname(p[, "Estimate"]),
                   unname(fitted(fit2)[1:5, "Estimate"]))
})

test_that("default_prior() does not read the response's values", {
  # priors:74: brms answers for a Beta() model on an rnorm() response
  set.seed(4)
  dd <- data.frame(y = rnorm(10), z = rnorm(10), g = rep(1:2, 5))
  gp <- default_prior(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()),
                      data = dd)
  expect_true(any(gp$class == "b" & gp$coef == "z" & gp$dpar == "phi"))
  # the fit still refuses the response, which is the absent case
  expect_error(frm(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()),
                   data = dd), "strictly in (0, 1)", fixed = TRUE)
})

test_that("se() on a family that cannot read it is refused before data", {
  # brm:81: brms refuses from the formula; frmtmb named the missing
  # column first
  set.seed(5)
  dd <- data.frame(y = rexp(10), x = rnorm(10))
  expect_error(frm(y | se(sei) ~ x, data = dd, family = weibull()),
               "'weibull' does not declare that it does", fixed = TRUE)
  # with the column there, the same refusal, from the frame as before
  dd$sei <- 1
  expect_error(frm(y | se(sei) ~ x, data = dd, family = weibull()),
               "'weibull' does not declare that it does", fixed = TRUE)
  dd$sei <- NULL
  # the absent case: a family that reads se() gets the missing column,
  # named by the addition-term rule (a variable must be in the data)
  expect_error(frm(y | se(sei) ~ x, data = dd, family = gaussian()),
               "reads `sei`, which is not a column of `data`", fixed = TRUE)
})

test_that("cs() inside a group-level term is refused by name", {
  # brm:116: R's "could not find function \"cs\"" before
  set.seed(6)
  dd <- data.frame(y = sample(1:3, 60, TRUE), x = rnorm(60),
                   g = factor(rep(1:6, 10)))
  dd$yf <- factor(dd$y)
  expect_error(frm(yf ~ x + (cs(x) | g), data = dd, family = categorical()),
               "cs() needs an sratio, cratio, or acat family", fixed = TRUE)
  expect_error(frm(y ~ x + (cs(x) | g), data = dd, family = sratio()),
               "is not supported: cs() is a population-level term here",
               fixed = TRUE)
})

test_that("family(resp = ) and print(links = ) are brms's", {
  # brmsfit-methods:283, :284, :285
  set.seed(7)
  dd <- data.frame(y1 = rnorm(40), y2 = rpois(40, 3), x = rnorm(40))
  fit <- frm(bf(y1 ~ x) + gaussian() + bf(y2 ~ x) + poisson(), data = dd)
  expect_identical(family(fit, resp = "y2")[["family"]], "poisson")
  expect_named(family(fit, resp = c("y1", "y2")), c("y1", "y2"))
  expect_error(family(fit, resp = "y3"), "Unknown response: 'y3'",
               fixed = TRUE)
  out <- paste(utils::capture.output(print(student(), links = TRUE)),
               collapse = "\n")
  expect_match(out, "Link function of 'nu' (if predicted): logm1",
               fixed = TRUE)
  out2 <- paste(utils::capture.output(print(student(), links = "sigma")),
                collapse = "\n")
  expect_false(grepl("'nu'", out2, fixed = TRUE))
  # the absent case: no link lines unless asked
  out0 <- paste(utils::capture.output(print(student())), collapse = "\n")
  expect_false(grepl("if predicted", out0, fixed = TRUE))
  mx <- paste(utils::capture.output(print(mixture(gaussian(),
                                                  exponential()))),
              collapse = "\n")
  expect_match(mx, "Mixture.*gaussian.*exponential")
  expect_error(print(student(), links = 3), "`links` must be")
})

test_that("mixture(order = ) is validated as brms validates it", {
  # families:102
  expect_error(mixture(poisson, binomial, order = "x"),
               "Argument 'order' is invalid", fixed = TRUE)
  expect_error(mixture(gaussian(), gaussian(), order = "mu"),
               "is not supported", fixed = TRUE)
  expect_error(mixture(gaussian(), gaussian(), order = TRUE),
               "is not supported", fixed = TRUE)
  expect_identical(mixture(gaussian(), gaussian(), order = "none"),
                   mixture(gaussian(), gaussian()))
  expect_identical(mixture(gaussian(), gaussian(), order = FALSE),
                   mixture(gaussian(), gaussian()))
})

test_that("ranef() takes brms's pars and groups", {
  # brmsfit-methods:814, :817
  set.seed(8)
  n <- 120
  dd <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)),
                   h = factor(rep(1:6, 20)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(12)[dd$g] + rnorm(6)[dd$h])
  fit <- frm(y ~ x + (1 + x | g) + (1 | h), data = dd)
  all <- ranef(fit)
  r <- ranef(fit, pars = "x")
  expect_named(r, "g")
  expect_identical(colnames(r$g), "x")
  expect_identical(unname(r$g[, "x"]), unname(all$g[, "x"]))
  expect_named(ranef(fit, groups = "h"), "h")
  expect_length(ranef(fit, groups = "a"), 0L)
  expect_error(ranef(fit, pars = 1), "`pars` must be a character vector")
})

test_that("fitted() takes sample_new_levels = \"gaussian\"", {
  # brmsfit-methods:317
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, dd$x + rnorm(6)[dd$g])
  fit <- frm(y ~ x + (1 | g), data = dd)
  nd <- data.frame(x = 0, g = "new")
  expect_identical(fitted(fit, newdata = nd, allow_new_levels = TRUE,
                          sample_new_levels = "gaussian"),
                   fitted(fit, newdata = nd, allow_new_levels = TRUE))
  expect_error(fitted(fit, newdata = nd, allow_new_levels = TRUE,
                      sample_new_levels = "old_levels"),
               "honors sample_new_levels = \"gaussian\" only", fixed = TRUE)
})

test_that("residuals() takes newdata and answers a multivariate fit", {
  # brmsfit-methods:828, :832, :840, :841
  set.seed(10)
  n <- 60
  dd <- data.frame(x = rnorm(n), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(n, 1 + dd$x + rnorm(6)[dd$g])
  fit <- frm(y ~ x + (1 | g), data = dd)
  nd <- dd[1:10, ]
  r <- residuals(fit, newdata = nd)
  f <- fitted(fit, newdata = nd)
  expect_identical(dim(r), c(10L, 4L))
  expect_identical(r[, "Estimate"], nd$y - f[, "Estimate"])
  expect_identical(r[, "Est.Error"], f[, "Est.Error"])
  # on the fitted rows it is the in-sample residual, arithmetic apart
  r0 <- residuals(fit)[1:10, ]
  expect_lt(max(abs(r - r0)), tol_eps(max(abs(r0))))
  nd$g <- factor(paste0("n", 1:10))
  expect_error(residuals(fit, newdata = nd), "New levels")
  expect_identical(dim(residuals(fit, newdata = nd,
                                 allow_new_levels = TRUE)), c(10L, 4L))
  expect_error(residuals(fit, newdata = nd[, "x", drop = FALSE],
                         allow_new_levels = TRUE),
               "needs the observed response")
  expect_error(residuals(fit, type = "pearson", newdata = dd[1:3, ]),
               "residuals(newdata = ) is defined", fixed = TRUE)

  dm <- data.frame(y1 = rnorm(n), y2 = rpois(n, 4), x = rnorm(n))
  fm <- frm(bf(y1 ~ x) + gaussian() + bf(y2 ~ x) + poisson(), data = dm)
  rm <- residuals(fm)
  expect_identical(dim(rm), c(as.integer(n), 4L, 2L))
  expect_identical(dimnames(rm)[[3]], c("y1", "y2"))
  expect_identical(rm[, , "y2"], residuals(fm, resp = "y2"))
  expect_identical(unname(rm[, "Estimate", "y1"]),
                   unname(dm$y1 - fitted(fm, resp = "y1")[, "Estimate"]))
  expect_error(residuals(fm, type = "osa"), "multivariate")
})

test_that("fitted(scale = \"linear\") on a cs() fit: a layer per threshold", {
  # brmsfit-methods:352. The shared predictor alone is the predictor at
  # no threshold of the model; brms returns eta + cs_k per threshold k
  set.seed(11)
  n <- 150
  dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y <- ifelse(runif(n) < plogis(-0.5 - 0.8 * dd$x1 - 0.3 * dd$x2), 1L, NA)
  p2 <- plogis(0.7 - 0.8 * dd$x1 + 0.4 * dd$x2)
  y[is.na(y)] <- ifelse(runif(sum(is.na(y))) < p2[is.na(y)], 2L, 3L)
  dd$y <- y
  fit <- frm(bf(y ~ x1 + cs(x2)), family = sratio(), data = dd)
  fl <- fitted(fit, scale = "linear")
  expect_identical(dim(fl), c(as.integer(n), 4L, 2L))
  # brms names the layers eta1, eta2 (dev/defects-brms-mvcs.R)
  expect_identical(dimnames(fl)[[3]], c("eta1", "eta2"))
  # IDENTITY: sratio's P(Y = 1) is F(tau_1 - eta_1) and P(Y = 2) is
  # (1 - P(Y = 1)) F(tau_2 - eta_2), with the density's own arithmetic
  fr <- fitted(fit)
  tau <- fixef(fit)[c("Intercept[1]", "Intercept[2]"), "Estimate"]
  p1 <- stats::plogis(tau[[1]] - fl[, "Estimate", 1])
  p2 <- (1 - p1) * stats::plogis(tau[[2]] - fl[, "Estimate", 2])
  expect_lt(max(abs(p1 - fr[, "Estimate", 1])), tol_eps(1))
  expect_lt(max(abs(p2 - fr[, "Estimate", 2])), tol_eps(1))
  expect_identical(dim(fitted(fit, newdata = dd[1, ], scale = "linear")),
                   c(1L, 4L, 2L))
  # the absent case: without cs() the linear predictor is one column
  f0 <- frm(bf(y ~ x1), family = sratio(), data = dd)
  expect_identical(dim(fitted(f0, scale = "linear")), c(as.integer(n), 4L))
})

test_that("summary() reports a gp() term's hyperparameters, as brms does", {
  # brmsfit-methods:904: the section had a printer and nothing to print
  set.seed(12)
  n <- 40
  dd <- data.frame(x = seq(0, 10, length.out = n))
  dd$y <- sin(dd$x) + rnorm(n, 0, 0.3)
  fit <- frm(y ~ gp(x), data = dd)
  s <- summary(fit)
  expect_identical(rownames(s$gp), c("sdgp(gpx)", "lscale(gpx)"))
  expect_output(print(s), "Gaussian Process Hyperparameters")
  expect_output(print(s), "sdgp(gpx)", fixed = TRUE)
  # brms's lscale is on inputs scaled by their largest distance, here 10
  cv <- confint_varcorr(fit)
  rg <- cv$estimate[cv$term == "range(gp)"]
  expect_lt(abs(s$gp["lscale(gpx)", "Estimate"] * 10 - rg), tol_eps(rg))
})

test_that("update(newdata = ) records the data's name, as brms does", {
  # brmsfit-methods:924
  set.seed(13)
  dd <- data.frame(x = rnorm(30))
  dd$y <- rnorm(30, dd$x)
  fit <- frm(y ~ x, data = dd)
  expect_identical(attr(fit$data, "data_name"), "dd")
  new_data <- dd[1:20, ]
  up <- update(fit, newdata = new_data)
  expect_identical(attr(up$data, "data_name"), "new_data")
  expect_output(print(summary(up)), "Data: new_data (Number", fixed = TRUE)
  # the absent case: a frame passed by value has no spelling to record
  byval <- do.call(frm, list(formula = y ~ x, data = dd))
  expect_null(attr(byval$data, "data_name"))
  expect_output(print(summary(byval)), "Data:  (Number", fixed = TRUE)
})

test_that("the prior table spells a smooth's column as brms does", {
  # priors:55, :59
  set.seed(14)
  dd <- data.frame(y = rnorm(30), x = rnorm(30), z = rnorm(30),
                   g = rep(1:3, 10))
  gp <- default_prior(y ~ z + s(x) + (1 | g), data = dd)
  expect_identical(gp$coef[gp$class == "b"], c("", "sx_1", "z"))
  gp2 <- default_prior(bf(y ~ lp, lp ~ z + s(x) + (1 | g), nl = TRUE),
                       data = dd)
  expect_identical(gp2$coef[gp2$class == "b"],
                   c("", "Intercept", "sx_1", "z"))
  # either spelling lands on the same coefficient
  lb <- function(co) {
    f <- frm(y ~ z + s(x), data = dd,
             prior = set_prior("", class = "b", coef = co, lb = 5))
    fixef_by_dpar(f)$mu[["s(x).fx1"]]
  }
  a <- lb("sx_1")
  expect_gte(a, 5)
  expect_identical(a, lb("s(x).fx1"))
  vp <- validate_prior(set_prior("normal(0, 1)", class = "b",
                                 coef = "s(x).fx1"),
                       y ~ z + s(x), data = dd)
  expect_identical(vp$prior[vp$class == "b" & vp$coef == "sx_1"],
                   "normal(0, 1)")
  expect_false(any(vp$coef == "s(x).fx1"))
})

test_that("predict() takes brms's sample_new_levels = \"old_levels\"", {
  # brmsfit-methods:772. An unseen level borrows ONE seen level's effect,
  # chosen once per call, as brms's get_new_rsamples() chooses it
  set.seed(15)
  ng <- 5
  dd <- data.frame(g = factor(rep(seq_len(ng), each = 30)))
  u <- c(-6, -3, 0, 3, 6)
  dd$y <- rnorm(nrow(dd), 10 + u[dd$g], 1)
  fit <- frm(y ~ 1 + (1 | g), data = dd)
  nd <- data.frame(g = c("new", "new"))
  nsim <- 2000
  d <- predict(fit, newdata = nd, allow_new_levels = TRUE,
               sample_new_levels = "old_levels", propagate_error = FALSE,
               ndraws = nsim, summary = FALSE)
  # both rows are the one unseen level, so both read one seen level
  m <- colMeans(d)
  se <- apply(d, 2, stats::sd) / sqrt(nsim)
  cand <- fixef(fit)["Intercept", "Estimate"] + ranef(fit)$g[, "Intercept"]
  gap <- sort(abs(m[1] - cand))
  expect_lt(gap[1], 5 * se[1])
  expect_gt(gap[2], 5 * se[1])
  expect_lt(abs(m[1] - m[2]), 5 * max(se))
  # the spread is the observation noise alone, not the between-group
  # variance "gaussian" adds
  dg <- predict(fit, newdata = nd, allow_new_levels = TRUE,
                propagate_error = FALSE, ndraws = nsim, summary = FALSE)
  expect_lt(stats::sd(d[, 1]), 0.5 * stats::sd(dg[, 1]))
  expect_error(predict(fit, newdata = nd, allow_new_levels = TRUE,
                       sample_new_levels = "uncertainty", ndraws = 10),
               "\"uncertainty\" mixes the two", fixed = TRUE)
})

test_that("\"old_levels\" picks the seen level at random, per call", {
  # the pick is one draw per call, as brms's sample(); a pick that is
  # always level 1 fails here
  set.seed(16)
  dd <- data.frame(g = factor(rep(1:5, each = 30)))
  dd$y <- rnorm(nrow(dd), 10 + c(-6, -3, 0, 3, 6)[dd$g], 1)
  fit <- frm(y ~ 1 + (1 | g), data = dd)
  cand <- fixef(fit)["Intercept", "Estimate"] + ranef(fit)$g[, "Intercept"]
  picked <- vapply(1:12, function(s) {
    set.seed(s)
    d <- predict(fit, newdata = data.frame(g = "new"),
                 allow_new_levels = TRUE, sample_new_levels = "old_levels",
                 propagate_error = FALSE, ndraws = 200, summary = FALSE)
    which.min(abs(mean(d) - cand))
  }, 1L)
  # 12 calls over 5 equally likely levels: one level every time has
  # probability 5 / 5^12
  expect_gt(length(unique(picked)), 1L)
})

test_that("a missing mi() response has no residual", {
  # punch round 1, B1: the frame holds the placeholder 0 at a missing
  # row, and residuals() reported 0 - mu there; brms's is NA
  set.seed(21)
  n <- 100
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- 1 + d$x + rnorm(n)
  d$ymi <- ifelse(runif(n) < 0.2, NA, d$y)
  miss <- is.na(d$ymi)
  fu <- frm(ymi | mi() ~ x, data = d)
  for (ty in c("response", "pearson", "deviance")) {
    rr <- residuals(fu, type = ty)
    r <- rr[, "Estimate"]
    expect_true(all(is.na(r[miss])), info = ty)
    expect_false(anyNA(r[!miss]), info = ty)
    # all four columns, as brms's (punch round 2, R1: Est.Error kept
    # the fitted value's standard error there)
    expect_true(all(is.na(rr[miss, ])), info = ty)
  }
  expect_false(anyNA(residuals(fu)[!miss, "Est.Error"]))
  # punch round 2, R2: TMB's own "'observation.name' must be in data
  # component" before
  expect_error(residuals(fu, type = "osa"),
               "not available for a response with mi()", fixed = TRUE)
  r <- residuals(fu)[, "Estimate"]
  expect_identical(unname(r[!miss]),
                   unname(d$ymi[!miss] - fitted(fu)[!miss, "Estimate"]))
  # the multivariate form, which residuals() refused before this lane
  d$xmi <- ifelse(runif(n) < 0.2, NA, d$x)
  mx <- is.na(d$xmi)
  fm <- frm(bf(y ~ mi(xmi)) + bf(xmi | mi() ~ z) + set_rescor(FALSE),
            data = d)
  rm <- residuals(fm)
  expect_true(all(is.na(rm[mx, "Estimate", "xmi"])))
  expect_false(anyNA(rm[!mx, "Estimate", "xmi"]))
  expect_false(anyNA(rm[, "Estimate", "y"]))
  # pp_check() drops the missing rows with brms's warning, where it
  # plotted the placeholders as observations
  skip_if_not_installed("bayesplot")
  allow_warnings(p <- pp_check(fu, ndraws = 5),
                 "NA responses are not shown",
                 require = "NA responses are not shown")
  expect_s3_class(p, "ggplot")
  skip_if_not_installed("DHARMa")
  dh <- dharma_residuals(fu, nsim = 20, seed = 1)
  expect_identical(length(dh$observedResponse), sum(!miss))
  expect_false(any(dh$observedResponse == 0))
})

test_that("fitted(scale = \"linear\") stacks a cs() response's layers", {
  # punch round 1, B2: a multivariate fit with an ordinal cs() response
  # died in array() on "length of 'dimnames'". brms 2.23.0 answers
  # 150 x 4 x 3 named yg, eta1, eta2 (dev/defects-brms-mvcs.R)
  set.seed(5)
  n <- 150
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$yo <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  d$yg <- d$z + rnorm(n)
  f <- frm(bf(yg ~ x) + bf(yo ~ x + cs(z), family = sratio()) +
             set_rescor(FALSE), data = d)
  a <- fitted(f, scale = "linear")
  expect_identical(dim(a), c(as.integer(n), 4L, 3L))
  expect_identical(dimnames(a)[[3]], c("yg", "eta1", "eta2"))
  expect_identical(a[, , "yg"], fitted(f, scale = "linear", resp = "yg"))
  yo <- fitted(f, scale = "linear", resp = "yo")
  expect_identical(a[, , "eta2"], yo[, , "eta2"])
})

test_that("the prior table refuses an ordinal response it cannot count", {
  # punch round 1, minor 3: the threshold rows are counted from the
  # response, so brms refuses, "Could not extract the number of
  # thresholds"; validate_prior() answers what default_prior() answers
  set.seed(17)
  dd <- data.frame(y = rnorm(10), z = rnorm(10))
  expect_error(default_prior(y ~ z, data = dd, family = cumulative()),
               "positive integers or ordered factors")
  dd$yo <- rep(1:3, length.out = 10)
  expect_true(any(default_prior(yo ~ z, data = dd,
                                family = cumulative())$class ==
                    "Intercept"))
  vp <- validate_prior(set_prior("normal(0, 1)", class = "b"),
                       bf(y ~ z, family = Beta()), data = dd)
  expect_identical(unique(vp$prior[vp$class == "b"]), "normal(0, 1)")
})
