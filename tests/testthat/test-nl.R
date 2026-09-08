test_that("nonlinear fixed-effects model matches nls", {
  set.seed(111)
  n <- 200
  x <- runif(n, 0, 5)
  y <- 2.5 * exp(-0.7 * x) + rnorm(n, 0, 0.15)
  dd <- data.frame(y = y, x = x)

  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(1, 0.3)))
  ref <- nls(y ~ a * exp(-b * x), data = dd, start = list(a = 1, b = 0.3))

  expect_lt(abs(fixef(fit)$a[[1]] - coef(ref)[["a"]]), 1e-4)
  expect_lt(abs(fixef(fit)$b[[1]] - coef(ref)[["b"]]), 1e-4)
  # ML sigma^2 = RSS/n at the same coefficients
  sig_ml <- sqrt(sum(residuals(ref)^2) / n)
  expect_lt(abs(exp(fixef(fit)$sigma[[1]]) - sig_ml), 1e-4)
})

test_that("nonlinear mixed model matches a hand-rolled reference", {
  set.seed(112)
  n_g <- 25; n_per <- 20
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- runif(n_g * n_per, 0, 5)
  a_g <- 2.5 + rnorm(n_g, 0, 0.5)
  y <- a_g[g] * exp(-0.7 * x) + rnorm(length(x), 0, 0.15)
  dd <- data.frame(y = y, x = x, g = g)

  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
               gaussian(),
             data = dd, start = list(beta = c(2, 0.5)))

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  nll_ref <- function(p) {
    nll <- -sum(RTMB::dnorm(p$u, 0, exp(p$lsd), log = TRUE))
    a <- p$a0 + p$u[gi]
    mu <- a * exp(-p$b0 * xv)
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(a0 = 2, b0 = 0.5, ls = 0, lsd = 0,
                              u = numeric(n_g)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  # nlpar random effects show up in ranef and VarCorr
  expect_length(VarCorr(fit), 1)
  expect_identical(dim(ranef(fit)[[1]]), c(as.integer(n_g), 1L))
})

test_that("nl prediction and post-processing", {
  set.seed(113)
  n <- 150
  x <- runif(n, 0, 5)
  dd <- data.frame(y = 2 * exp(-0.5 * x) + rnorm(n, 0, 0.1), x = x)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(1, 0.3)))

  expect_equal(predict(fit, newdata = dd), predict(fit), tolerance = 1e-8)
  expect_equal(fitted(fit), predict(fit, type = "response"),
               tolerance = 1e-8)
  a_hat <- predict(fit, dpar = "a")
  expect_lt(stats::sd(a_hat), 1e-10)   # intercept-only nlpar is constant
  nd <- data.frame(x = c(0, 1, 2))
  p <- predict(fit, newdata = nd)
  expect_equal(p[1], fixef(fit)$a[[1]], tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_error(predict(fit, se.fit = TRUE), "se.fit is not supported")
  expect_length(residuals(fit), n)
})

test_that("nl validation errors are clear", {
  expect_error(bf(y ~ a * exp(-b * x), nl = TRUE), "parameter formula")
  expect_error(frm(bf(y ~ a * exp(-b * x), a ~ 1, cc ~ 1, nl = TRUE) +
                     gaussian(),
                   data = NULL, dry_run = "spec"),
               "not used in the model formula")
})

# --- reserved nonlinear parameter names -------------------------------
# Four usability defects of hierarchical nl models, wt-api. A nonlinear
# parameter named after one of the family's own dpars used to die inside
# model.frame() with "object 'mu' not found", naming nothing; one named
# after a par-template component fits, but `start` then means the
# component and reported ITS length.

nl_reserved_data <- function(n_id = 8, seed = 7) {
  set.seed(seed)
  w <- seq(1, 6, by = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    data.frame(id = i, group = i %% 2, w = w, logw = log(w),
               I = rexp(length(w), rate = 1 / exp(1.2 - 1.5 * log(w))))
  }))
  d$id <- factor(d$id)
  d
}

test_that("a nonlinear parameter named after a family dpar is refused by name", {
  d <- nl_reserved_data()
  # the reported spelling: `mu ~ 1 + (1 | id)` alongside nl = TRUE
  expect_error(
    frm(bf(I ~ mu - chi * logw, mu ~ 1 + (1 | id), chi ~ 1 + group,
           nl = TRUE),
        family = exponential(link = "log"), data = d),
    "distributional parameter of family 'exponential'")
  expect_error(
    frm(bf(I ~ mu - chi * logw, mu ~ 1 + (1 | id), chi ~ 1 + group,
           nl = TRUE),
        family = exponential(link = "log"), data = d),
    "This family reserves: mu")
  # and at SPEC time, so par_template() refuses it before any fit
  expect_error(
    par_template(bf(I ~ mu - chi * logw, mu ~ 1 + (1 | id),
                    chi ~ 1 + group, nl = TRUE),
                 data = d, family = exponential(link = "log")),
    "cannot also be a nonlinear parameter")
  # any family, not just this one
  expect_error(
    frm(bf(I ~ mu * chi, mu ~ 1, chi ~ 1, nl = TRUE), gaussian(), data = d),
    "distributional parameter of family 'gaussian'")
})

test_that("a body that names its own parameter is refused, with the data checked first", {
  d <- nl_reserved_data()
  # no formula for mu and no column called mu: the old message was R's
  # own "object 'mu' not found" from eval(predvars, data, env)
  expect_error(
    frm(bf(I ~ mu - chi * logw, chi ~ 1 + group, nl = TRUE),
        family = exponential(link = "log"), data = d),
    "refers to 'mu' itself")
  # a REAL column of that name still wins, as it does for a dpar
  # reference, so this model keeps fitting
  d2 <- d
  d2$mu <- 1
  f <- frm(bf(I ~ mu * apo - chi * logw, apo ~ 1, chi ~ 1 + group,
              nl = TRUE),
           family = exponential(link = "log"), data = d2)
  # mu is computed by the body, so it contributes no coefficient block;
  # what matters is that the fit happened at all
  expect_setequal(names(fixef(f)), c("apo", "chi"))
  expect_true(all(is.finite(fixef(f, flatten = TRUE))))
  # an nlf() body that names ITSELF is the same fault under another
  # spelling
  expect_error(
    frm(bf(I ~ apo - chi * logw, apo ~ 1, chi ~ 1 + group, nl = TRUE) +
          nlf(sigma ~ sigma + 1),
        family = gaussian(), data = d),
    "refers to 'sigma' itself")
})

test_that("a body reading ANOTHER dpar's value is untouched", {
  # the deliberate variance-function extension: `sigma` in mu's body is
  # that parameter's per-row value, not a nonlinear parameter, and the
  # refusals above must not reach it
  set.seed(2)
  dd <- data.frame(x = rnorm(120))
  dd$y <- 2 + dd$x + rnorm(120)
  f <- frm(bf(y ~ sigma * x + a, a ~ 1, nl = TRUE), gaussian(), data = dd)
  expect_true(is.finite(fixef(f)$a[[1]]))
})

test_that("a nonlinear parameter named after a template component fits, and start names the collision", {
  d <- nl_reserved_data()
  # one test per par-template component that a plain hierarchical model
  # carries. Each name FITS - the refusal is only about `start`.
  for (nm in c("beta", "b", "theta")) {
    body <- stats::as.formula(paste0("I ~ ", nm, " - chi * logw"))
    par <- stats::as.formula(paste0(nm, " ~ 1 + (1 | id)"))
    f <- frm(bf(body, par, chi ~ 1 + group, nl = TRUE),
             family = exponential(link = "log"), data = d)
    expect_true(paste0(nm, "_(Intercept)") %in% names(fixef(f, flatten = TRUE)),
                info = nm)
    st <- list(1)
    names(st) <- nm
    # `start$<nm>` sets the COMPONENT. Whether that is a length error or
    # a silent success depends on the component's length, so both paths
    # have to carry the explanation.
    seen <- character(0)
    msg <- tryCatch({
      withCallingHandlers(
        frm(bf(body, par, chi ~ 1 + group, nl = TRUE),
            family = exponential(link = "log"), data = d, start = st),
        warning = function(w) {
          seen <<- c(seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      paste(seen, collapse = " ")
    }, error = function(e) conditionMessage(e))
    expect_match(msg, "not the nonlinear parameter", info = nm)
    expect_match(msg, "par_template\\(\\) lists both", info = nm)
  }
})

test_that("every par-template component name is spelled out by the collision message", {
  # thetaac, thetar and miss need an autocorrelation term, a residual
  # correlation and an imputed column to appear in a template at all;
  # the message that names them is unit-tested instead of fitting three
  # more models for one string each
  tpl <- list(beta = c(`z_(Intercept)` = 0), betad = numeric(0),
              b = numeric(0), theta = numeric(0), thetaac = numeric(0),
              thetar = numeric(0), miss = numeric(0))
  expect_match(frmtmb:::nl_start_collision_msg("beta", tpl),
               "fixed-effect coefficients")
  expect_match(frmtmb:::nl_start_collision_msg("betad", tpl),
               "distributional parameters")
  expect_match(frmtmb:::nl_start_collision_msg("b", tpl),
               "random-effect vector")
  expect_match(frmtmb:::nl_start_collision_msg("theta", tpl),
               "covariance parameters")
  expect_match(frmtmb:::nl_start_collision_msg("thetaac", tpl),
               "autocorrelation parameters")
  expect_match(frmtmb:::nl_start_collision_msg("thetar", tpl),
               "residual-correlation parameters")
  expect_match(frmtmb:::nl_start_collision_msg("miss", tpl),
               "imputed missing values")
})

test_that("newparams carries the same collision message, in its own spelling", {
  d <- nl_reserved_data()
  msg <- tryCatch(
    frm_simulate(bf(I ~ b - chi * logw, b ~ 1 + (1 | id), chi ~ 1 + group,
                    nl = TRUE),
                 data = d, family = exponential(link = "log"),
                 newparams = list(beta = c(1, 1.5, 0.2), b = 1,
                                  theta = 0.3),
                 nsim = 1, seed = 1),
    error = function(e) conditionMessage(e))
  expect_match(msg, "not the nonlinear parameter")
  # the message answers the argument the caller actually used
  expect_match(msg, "`newparams[$]b`")
  expect_false(grepl("start$", msg, fixed = TRUE))
})

test_that("the collision message keeps saying start$ for start", {
  d <- nl_reserved_data()
  msg <- tryCatch(
    frm(bf(I ~ b - chi * logw, b ~ 1 + (1 | id), chi ~ 1 + group,
           nl = TRUE),
        family = exponential(link = "log"), data = d, start = list(b = 1)),
    error = function(e) conditionMessage(e))
  expect_match(msg, "`start[$]b`")
  expect_false(grepl("newparams", msg, fixed = TRUE))
})
