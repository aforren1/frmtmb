# A term that calls a function with a covariance-structure name, such as
# exp(x), is protected from reformulas::splitForm() by an alias
# (`.frm_exp`, dev/rtmb-pitfalls.md item 9). Through 0.74.0 the alias
# reached model.matrix(), so the coefficient of exp(x) was `.frm_expx`
# in fixef(), vcov(), summary(), coef() and VarCorr(), `b_.frm_expx` in
# variables() and the draws, and `.frm_exp(x)` in confint() and
# default_prior(); brms 2.23.0 names it `expx` and `b_expx`, so
# set_prior(coef = "expx") and hypothesis(fit, "expx > 0") were refused.
# The prior table also listed every transformed column under its
# model.matrix() name (`log(xp)`, `I(x^2)`, `poly(x, 2)1`), where brms
# lists its rename() (`logxp`, `IxE2`, `polyx21`), and refused brms's
# spelling for any column whose rename is more than dropped parentheses.
#
# brms's names are computed here from brms itself (standata() and
# default_prior(), no compile). The table of every term against brms,
# with each interface, is dev/brmsnames-1024-findings.md
# (dev/brmsnames-1024-table.R). Seen to fail on 0.74.0 (rellib-r12).

bns_data <- function() {
  set.seed(1024)
  n <- 120
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  xp = stats::runif(n, 0.5, 2),
                  k = sample(1:4, n, replace = TRUE),
                  f = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
                  g = factor(rep(1:12, length.out = n)),
                  e = stats::runif(n, 1, 3))
  d$y <- 1 + 0.4 * exp(d$x / 2) + 0.3 * d$z + 0.2 * d$k +
    stats::rnorm(12, 0, 0.3)[d$g] + stats::rnorm(n, 0, exp(0.1 * d$z))
  d$y2 <- 0.5 * d$x + stats::rnorm(n)
  d$yc <- stats::rpois(n, d$e * exp(0.2 + 0.3 * d$x))
  d$yn <- 2 * exp(0.3 * d$x) + stats::rnorm(n, 0, 0.2)
  d
}

# Each case: the formula text, read by brms's bf() and by frmtmb's, and
# the family. log(), sqrt(), I(), poly(), scale() and factor() were
# already named as brms names them in fixef(); they are here as the
# case where only the prior table changes.
bns_cases <- list(
  list("bf(y ~ exp(x))", "gaussian"),
  list("bf(y ~ exp(x):z)", "gaussian"),
  list("bf(y ~ exp(x) * f)", "gaussian"),
  list("bf(y ~ I(exp(x)))", "gaussian"),
  list("bf(y ~ exp(x / 2))", "gaussian"),
  list("bf(y ~ exp(x) + log(xp) + sqrt(xp))", "gaussian"),
  list("bf(y ~ I(x^2) + scale(z))", "gaussian"),
  list("bf(y ~ poly(z, 2, raw = TRUE) + factor(k))", "gaussian"),
  list("bf(y ~ x, sigma ~ exp(z))", "gaussian"),
  list("bf(y ~ a + b * x, a ~ exp(z), b ~ 1, nl = TRUE)", "gaussian"),
  list("bf(y ~ exp(x) + (1 + exp(x) | g))", "gaussian"),
  list("bf(y ~ exp(x) + (1 + exp(x) || g))", "gaussian"),
  list("mvbf(bf(y ~ exp(x)), bf(y2 ~ log(xp)), rescor = FALSE)", "gaussian"),
  list("bf(yc ~ exp(x) + offset(log(e)))", "poisson"),
  # brms's prefix is dpar, resp, nlpar: yn_a_expz, sigma_yn_z
  list(paste("mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + exp(z), b ~ 1,",
             "sigma ~ z, nl = TRUE), bf(y ~ x + exp(z)), rescor = FALSE)"),
       "gaussian")
)

bns_form <- function(txt, pkg) {
  env <- list2env(list(bf = getExportedValue(pkg, "bf"),
                       mvbf = getExportedValue(pkg, "mvbf")),
                  parent = globalenv())
  eval(str2lang(txt), env)
}

# brms's fixef() names: each design matrix of standata(), its columns
# behind the predictor's prefix (X_sigma is sigma_, X_a is a_, X_y2 is
# y2_, X_yn_a is yn_a_), which is how brms's rename_pars() names b
bns_brms_fixef <- function(txt, fam, d) {
  sd <- brms_standata(bns_form(txt, "brms"), data = d, family = fam)
  xs <- grep("^X(_.+)?$", names(sd), value = TRUE)
  unlist(lapply(xs, function(nm) {
    pre <- sub("^X_?", "", nm)
    cn <- colnames(sd[[nm]])
    if (nzchar(pre)) paste0(pre, "_", cn) else cn
  }), use.names = FALSE)
}

bns_prior_rows <- function(p) {
  p <- as.data.frame(p)
  p <- p[nzchar(p$coef) & p$class %in% c("b", "sd"), , drop = FALSE]
  p[, c("class", "coef", "group", "resp", "dpar", "nlpar")]
}

bns_key <- function(p) do.call(paste, c(unname(as.list(p)), sep = "|"))

bns_fit <- function(txt, fam, d) {
  allow_boundary(frm(bns_form(txt, "frmtmb"), data = d, family = fam))
}

test_that("population-level names are brms's in every accessor", {
  skip_unless_brms()
  d <- bns_data()
  for (cs in bns_cases) {
    want <- bns_brms_fixef(cs[[1]], cs[[2]], d)
    fit <- bns_fit(cs[[1]], cs[[2]], d)
    lab <- cs[[1]]
    fx <- rownames(fixef(fit))
    expect_setequal(fx, want)
    expect_identical(rownames(vcov(fit)), fx, label = lab)
    expect_identical(rownames(summary(fit)$fixed), fx, label = lab)
    expect_true(all(paste0("b_", want) %in% variables(fit)), label = lab)
    for (nm in want) {
      h <- hypothesis(fit, paste0(nm, " = 0"))
      expect_s3_class(h, "frmtmb_hypothesis")
    }
  }
})

test_that("the prior table lists brms's coefficients and takes them", {
  skip_unless_brms()
  d <- bns_data()
  for (cs in bns_cases) {
    lab <- cs[[1]]
    bp <- bns_prior_rows(suppressMessages(brms::default_prior(
      bns_form(cs[[1]], "brms"), data = d, family = cs[[2]])))
    fp <- bns_prior_rows(default_prior(bns_form(cs[[1]], "frmtmb"),
                                       data = d, family = cs[[2]]))
    expect_setequal(bns_key(fp), bns_key(bp))
    for (i in seq_len(nrow(bp))) {
      r <- bp[i, ]
      a <- list(prior = if (r$class == "sd") "student_t(3, 0, 5)" else
        "normal(0, 5)", class = r$class, coef = r$coef)
      for (nm in c("group", "resp", "dpar", "nlpar")) {
        if (nzchar(r[[nm]])) a[[nm]] <- r[[nm]]
      }
      v <- as.data.frame(validate_prior(do.call(set_prior, a),
                                        bns_form(cs[[1]], "frmtmb"),
                                        data = d, family = cs[[2]]))
      expect_true(any(v$source == "user" & v$coef == r$coef &
                        v$class == r$class),
                  label = paste(lab, r$class, r$coef))
    }
  }
})

test_that("group-level names of exp(x) are brms's", {
  skip_unless_brms()
  d <- bns_data()
  for (txt in c("bf(y ~ exp(x) + (1 + exp(x) | g))",
                "bf(y ~ exp(x) + (1 + exp(x) || g))")) {
    bp <- as.data.frame(suppressMessages(brms::default_prior(
      bns_form(txt, "brms"), data = d, family = "gaussian")))
    sdc <- bp$coef[bp$class == "sd" & nzchar(bp$coef)]
    fit <- bns_fit(txt, "gaussian", d)
    expect_true(all(paste0("sd_g__", sdc) %in% variables(fit)), label = txt)
    expect_setequal(rownames(VarCorr(fit)$g$sd), sdc)
    expect_setequal(colnames(coef(fit)$g), c("Intercept", "expx"))
    rf <- ranef(fit)$g
    expect_setequal(if (length(dim(rf)) == 3L) dimnames(rf)[[3]] else
      colnames(rf), sdc)
    # the slope's sd can sit at its boundary on these data; the name is
    # what is tested, and hypothesis() says so when its SE is NaN
    h <- allow_warnings(hypothesis(fit, "g__expx = 0", class = "sd"),
                        "a direction the fit does not determine")
    expect_s3_class(h, "frmtmb_hypothesis")
  }
})

test_that("no placeholder reaches a name a user sees", {
  d <- bns_data()
  for (cs in bns_cases) {
    fit <- bns_fit(cs[[1]], cs[[2]], d)
    seen <- c(rownames(fixef(fit)), variables(fit), rownames(confint(fit)),
              names(unlist(lapply(coef(fit), function(x) {
                if (is.data.frame(x)) colnames(x) else names(x)
              }))),
              as.data.frame(default_prior(fit))$coef,
              names(fixef(fit, flatten = TRUE)))
    expect_false(any(grepl(".frm_", seen, fixed = TRUE)), label = cs[[1]])
  }
})

# The design column keeps R's spelling where frmtmb uses it internally
# (confint() rows, fixef(flatten = TRUE)), and that spelling is no
# longer the alias.
test_that("confint() names the design column as written", {
  d <- bns_data()
  fit <- bns_fit("bf(y ~ exp(x) * f)", "gaussian", d)
  expect_identical(rownames(confint(fit))[1:6],
                   c("(Intercept)", "exp(x)", "fb", "fc", "exp(x):fb",
                     "exp(x):fc"))
  # brms's spelling and the bare one address the same row
  expect_identical(rownames(confint(fit, parm = "b_expx")), "exp(x)")
  expect_identical(rownames(confint(fit, parm = "expx:fb")), "exp(x):fb")
})

# brms centers a predictor's intercept and stores it ahead of the b
# vector, but a nonlinear parameter's sub-formula is not centered: its
# intercept is an element of its b vector. brms 2.23.0's fixef() order
# for this model, from a brmsfit (dev/brmsnames-1024-brms.R, case 39,
# data of dev/brmsnames-1024-cases.R): through 0.74.0 frmtmb put
# a_Intercept and b_Intercept ahead of everything.
test_that("a nonlinear parameter's intercept keeps its place, as in brms", {
  d <- bns_data()
  fit <- bns_fit(paste("bf(y ~ a + b * x, a ~ exp(z), b ~ 1 + z,",
                       "sigma ~ z, nl = TRUE)"), "gaussian", d)
  want <- c("sigma_Intercept", "a_Intercept", "a_expz", "b_Intercept",
            "b_z", "sigma_z")
  expect_identical(rownames(fixef(fit)), want)
  expect_identical(rownames(vcov(fit)), want)
  expect_identical(rownames(summary(fit)$fixed), want)
})

# A multivariate model's nonlinear parameter: brms's prefix order is
# dpar, resp, nlpar (vars_prefix()), so yn_a_Intercept and
# sd_g__yn_a_Intercept, where frmtmb through 0.74.0 wrote a_yn_Intercept
# and refused brms's names in hypothesis(). The fixef() order is brms
# 2.23.0's from a brmsfit (dev/brmsnames-1024-brms.R, case 40): response
# by response, the centered intercepts first.
test_that("a multivariate nonlinear parameter is named as brms names it", {
  d <- bns_data()
  fit <- bns_fit(paste("mvbf(bf(yn ~ a * exp(b * x),",
                       "a ~ 1 + exp(z) + (1 | g), b ~ 1, sigma ~ z,",
                       "nl = TRUE), bf(y ~ x + exp(z)), rescor = FALSE)"),
                 "gaussian", d)
  want <- c("sigma_yn_Intercept", "y_Intercept", "yn_a_Intercept",
            "yn_a_expz", "yn_b_Intercept", "sigma_yn_z", "y_x", "y_expz")
  expect_identical(rownames(fixef(fit)), want)
  expect_true(all(c(paste0("b_", want), "sd_g__yn_a_Intercept") %in%
                    variables(fit)))
  for (nm in c("yn_a_Intercept", "yn_a_expz", "yn_b_Intercept")) {
    expect_s3_class(hypothesis(fit, paste0(nm, " = 0")), "frmtmb_hypothesis")
  }
  h <- allow_warnings(hypothesis(fit, "g__yn_a_Intercept = 0", class = "sd"),
                      "a direction the fit does not determine")
  expect_s3_class(h, "frmtmb_hypothesis")
  expect_identical(rownames(VarCorr(fit)$g$sd), "yn_a_Intercept")
})

# Two answers the alias got wrong, beyond names: drop1() read the term
# as one with no columns (Df 0, the full model's AIC), and a re_formula
# that names the fit's own term matched nothing. Seen to fail on 0.74.0.
test_that("drop1() and step() see an exp() term", {
  set.seed(3)
  n <- 200
  d <- data.frame(x = stats::runif(n), z = stats::rnorm(n))
  d$y <- 1 + 0.8 * exp(d$x) + 0.05 * d$z + stats::rnorm(n)
  fit <- frm(y ~ exp(x) + z, data = d)
  d1 <- drop1(fit)
  expect_identical(rownames(d1), c("<none>", "exp(x)", "z"))
  expect_equal(d1["exp(x)", "Df"], 1)
  red <- AIC(frm(y ~ z, data = d))
  expect_lt(abs(d1["exp(x)", "AIC"] - red), 1e-6 * abs(red))
  # dropping z lowers the AIC, dropping exp(x) raises it
  s <- step(fit, trace = 0)
  expect_identical(rownames(fixef(s)), c("Intercept", "expx"))
})

test_that("re_formula can name a group-level exp() term", {
  set.seed(77)
  n <- 300
  d <- data.frame(x = stats::runif(n, 0.2, 2),
                  g = factor(rep(1:20, each = 15)))
  d$y <- 1 + 0.2 * exp(d$x) + stats::rnorm(20, 0, 0.4)[d$g] +
    stats::rnorm(20, 0, 0.5)[d$g] * exp(d$x) + stats::rnorm(n, 0, 0.3)
  fit <- allow_boundary(frm(y ~ exp(x) + (1 + exp(x) | g), data = d))
  nd <- d[c(1, 20, 40), ]
  full <- fitted(fit, newdata = nd)
  same <- fitted(fit, newdata = nd, re_formula = ~ (1 + exp(x) | g))
  expect_identical(same, full)
  icpt <- fitted(fit, newdata = nd, re_formula = ~ (1 | g))
  expect_false(isTRUE(all.equal(icpt, full)))
})

# The design column's spelling keeps working beside brms's, for every
# class that takes a coef, and the prior lands on the coefficient it
# names: a normal(3, 0.001) moves that one coefficient to 3.
bns_prior_data <- function() {
  set.seed(31)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  xp = stats::runif(n, 0.5, 2),
                  k = factor(sample(1:3, n, replace = TRUE)),
                  g = factor(rep(1:10, length.out = n)))
  d$y <- 1 + d$x + 0.3 * d$z^2 + 0.5 * log(d$xp) +
    stats::rnorm(10)[d$g] + stats::rnorm(n)
  d
}

test_that("set_prior(coef =) takes the design column and brms's name", {
  d <- bns_prior_data()
  form <- bf(y ~ poly(x, 2) + I(z^2) + log(xp) + factor(k))
  n_rows <- nrow(as.data.frame(default_prior(form, data = d)))
  pairs <- list(c("poly(x, 2)1", "polyx21"), c("I(z^2)", "IzE2"),
                c("log(xp)", "logxp"), c("factor(k)2", "factork2"))
  for (pr in pairs) {
    for (spell in pr) {
      p <- set_prior("normal(3, 0.001)", class = "b", coef = spell)
      v <- as.data.frame(validate_prior(p, form, data = d))
      expect_identical(nrow(v), n_rows)
      expect_identical(v$coef[v$source == "user"], pr[2], label = spell)
      fe <- fixef(frm(form, data = d, prior = p))[, "Estimate"]
      expect_identical(names(fe)[abs(fe - 3) < 0.01 * 3], pr[2],
                       label = spell)
      ps <- prior_string("normal(3, 0.001)", class = "b", coef = spell)
      fe2 <- fixef(frm(form, data = d, prior = ps))[, "Estimate"]
      expect_identical(fe2, fe)
    }
  }
  # the two spellings name one slot, so writing both is a duplicate
  expect_error(validate_prior(
    set_prior("normal(0, 1)", class = "b", coef = "I(z^2)") +
      set_prior("normal(0, 2)", class = "b", coef = "IzE2"),
    form, data = d), "Duplicated prior")
})

# A smooth's unpenalized column, s(x).fx1, and brms's sx_1 are one slot
# too. Through 0.74.0 they were two, and the later prior replaced the
# earlier one without a word (dev/brmsnames-1024-smoothdup.R).
test_that("a smooth column's two spellings are one prior slot", {
  d <- bns_prior_data()
  fs <- bf(y ~ s(x) + z)
  expect_error(validate_prior(
    set_prior("normal(3, 0.001)", class = "b", coef = "s(x).fx1") +
      set_prior("normal(-3, 0.001)", class = "b", coef = "sx_1"),
    fs, data = d), "Duplicated prior")
  expect_error(allow_boundary(frm(fs, data = d, prior =
    set_prior("normal(3, 0.001)", class = "b", coef = "s(x).fx1") +
      set_prior("normal(-3, 0.001)", class = "b", coef = "sx_1"))),
    "Duplicated prior")
})

test_that("class sd and cs() take the design column's spelling too", {
  d <- bns_prior_data()
  for (cs in list(c("y ~ I(x^2) + (1 + I(x^2) | g)", "I(x^2)", "IxE2"),
                  c("y ~ poly(x, 2) + (1 + poly(x, 2) | g)", "poly(x, 2)2",
                    "polyx22"))) {
    form <- bf(stats::as.formula(cs[1]))
    for (spell in cs[2:3]) {
      p <- set_prior("normal(2, 0.001)", class = "sd", coef = spell,
                     group = "g")
      v <- as.data.frame(validate_prior(p, form, data = d))
      expect_identical(v$coef[v$source == "user"], cs[3], label = spell)
      # a prior this tight on one sd can leave a correlation's gradient
      # short of the tolerance; which sd it moved is the point here
      fit <- allow_warnings(allow_boundary(frm(form, data = d, prior = p)),
                            c("Standard errors are not available",
                              "standard error",
                              "Large maximum absolute gradient"))
      sdv <- VarCorr(fit)$g$sd[, "Estimate"]
      expect_identical(names(sdv)[abs(sdv - 2) < 0.01 * 2], cs[3],
                       label = spell)
    }
  }
  d$yo <- cut(d$x + stats::rlogis(nrow(d)), c(-Inf, -1, 0, 1, Inf),
              labels = FALSE)
  form <- bf(yo ~ cs(I(x^2)) + z)
  for (spell in c("I(x^2)", "IxE2")) {
    v <- as.data.frame(validate_prior(
      set_prior("normal(0, 1)", class = "b", coef = spell), form, data = d,
      family = acat()))
    expect_identical(v$coef[v$source == "user"], "IxE2", label = spell)
  }
})

# exp() in the second response's nonlinear formula: through 0.74.0 the
# fit stopped with 'could not find function ".frm_exp"', because the
# alias was bound where only the first response's frame looked
# (dev/brmsnames-1024-mvnl2.R).
test_that("a second nonlinear response may call exp()", {
  d <- bns_data()
  fit <- bns_fit(paste("mvbf(bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1,",
                       "nl = TRUE), bf(y ~ c + d * x, c ~ exp(z), d ~ 1,",
                       "nl = TRUE), rescor = FALSE)"), "gaussian", d)
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_true("y_c_expz" %in% rownames(fixef(fit)))
})
