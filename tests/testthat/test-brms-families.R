# Family objects and links, held to brms 2.23.0 (lane wt-famlink).
#
# The first block is brms's own tests.families.R, the assertions that
# transfer, verbatim. dev/famlink-port-ledger.R runs EVERY assertion of
# that file against a library and dev/famlink-findings.md records, per
# assertion, pass, deliberate divergence or cannot transfer.

test_that("brms tests.families.R: family functions (ported verbatim)", {
  expect_equal(student(identity)$link, "identity")
  expect_equal(student()$link, "identity")
  expect_equal(bernoulli(logit)$link, "logit")
  expect_error(bernoulli("sqrt"), "bernoulli")
  expect_equal(negbinomial(sqrt)$link, "sqrt")
  expect_error(negbinomial(inverse), "inverse")
  expect_equal(geometric(identity)$link, "identity")
  expect_error(geometric("inv"), "geometric")
  expect_equal(exponential(log)$link, "log")
  expect_error(exponential("cloglog"), "exponential")
  expect_equal(weibull()$family, "weibull")
  expect_error(weibull(sqrt), "weibull")
  expect_equal(Beta("probit")$link, "probit")
  expect_error(Beta("1/mu^2"), "beta")
  expect_equal(hurdle_poisson()$link, "log")
  expect_equal(hurdle_gamma()$family, "hurdle_gamma")
  expect_error(hurdle_gamma(sqrt), "hurdle_gamma")
  expect_equal(zero_inflated_poisson(log)$link, "log")
  expect_error(zero_inflated_poisson(list(1)), "zero_inflated_poisson")
  expect_equal(zero_inflated_negbinomial("log")$link, "log")
  expect_error(zero_inflated_negbinomial("logit"),
               "zero_inflated_negbinomial")
  expect_equal(zero_inflated_beta(logit)$family, "zero_inflated_beta")
  expect_equal(zero_inflated_binomial()$link_zi, "logit",
               ignore_attr = TRUE)
  expect_error(zero_inflated_binomial(y ~ x), "zero_inflated_binomial")
  expect_equal(categorical()$link, "logit")
  expect_error(categorical(probit), "probit")
  expect_equal(cumulative(cauchit)$family, "cumulative")
  expect_equal(sratio(probit_approx)$link, "probit_approx")
  expect_equal(cratio("cloglog")$family, "cratio")
  expect_equal(brmsfamily("gaussian", inverse)$link, "inverse")
  expect_equal(brmsfamily("geometric", "identity")$family, "geometric")
  expect_equal(brmsfamily("zi_poisson")$link_zi, "logit")
  expect_error(weibull(link_shape = "logit"),
               "'logit' is not a supported link for parameter 'shape'")
  expect_equal(beta_binomial()$link, "logit")
  expect_equal(beta_binomial("probit")$link, "probit")
  expect_equal(beta_binomial()$link_phi, "log")
  expect_error(beta_binomial("log"))
  expect_error(beta_binomial(link_phi = "logit"))
})

test_that("brms tests.families.R: print and mixture (ported verbatim)", {
  expect_output(print(weibull()), "Family: weibull \nLink function: log")
  expect_error(mixture(gaussian, "x"), "x is not a supported family")
  expect_error(mixture(poisson(), categorical()),
               "Some of the families are not allowed in mixture models")
  expect_error(mixture(poisson, "cumulative"),
               "Cannot mix ordinal and non-ordinal families")
  expect_error(mixture(lognormal, exgaussian, poisson()),
               "Cannot mix families with real and integer support")
  expect_error(mixture(lognormal),
               "Expecting at least 2 mixture components")
})

test_that("brms tests.brm.R and tests.standata.R: families, links, responses", {
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
  # brm() is frm(); standata() is frm(dry_run = "frame")
  expect_error(frm(y ~ x, dat, family = poisson("inverse")),
               "'inverse' is not a supported link for family 'poisson'")
  expect_error(frm(y ~ x, dat, family = c("weibull", "sqrt")),
               "'sqrt' is not a supported link for family 'weibull'")
  expect_error(frm(y ~ x, dat, family = c("categorical", "probit")),
               "'probit' is not a supported link for family 'categorical'")
  expect_error(frm(y ~ x, dat, family = "ordinal"),
               "ordinal is not a supported family")
  expect_error(frm(y ~ 1, data = data.frame(y = factor(-1:1)),
                   family = "cratio", dry_run = "frame"),
               paste("Family 'cratio' requires either positive integers",
                     "or ordered factors"))
  expect_error(frm(y ~ 1, data = data.frame(y = rep(0.5:7.5), 2),
                   family = "sratio", dry_run = "frame"),
               paste("Family 'sratio' requires either positive integers",
                     "or ordered factors"))
})

# ------------------------------------------------ the mean-link sets

test_that("the committed mean-link sets are the ones brms has now", {
  skip_if_not_installed("brms")
  ns <- asNamespace("brms")
  src <- frmtmb:::brms_mu_link_source
  live <- lapply(src, function(b) {
    get(paste0(".family_", b), envir = ns)()$links
  })
  expect_identical(live, frmtmb:::brms_mu_links)
  # a guard that compared nothing would pass on an empty table
  expect_gt(length(live), 30L)
  src <- src[setdiff(names(src), frmtmb:::brms_mu_link_analogs)]
  barred <- names(src)[vapply(src, function(b) {
    isTRUE(brms:::no_mixture(brms::brmsfamily(b)))
  }, NA)]
  expect_setequal(barred, frmtmb:::brms_no_mixture)
  expect_true("zero_inflated_beta" %in% barred)
})

# The registry constructor for each keyed family name.
famlink_ctor <- function(fam) frmtmb:::family_ctor(fam)

test_that("every (family, link) pair brms accepts constructs", {
  sets <- frmtmb:::brms_mu_links
  built <- 0L
  for (fam in setdiff(names(sets), "multinomial")) {
    ctor <- famlink_ctor(fam)
    for (lk in sets[[fam]]) {
      if (fam == "acat" && lk != "logit") next  # deliberate, below
      f <- ctor(link = lk)
      expect_identical(f$link, lk, label = paste(fam, lk))
      built <- built + 1L
    }
  }
  expect_identical(built, sum(lengths(sets[setdiff(names(sets),
                                                   "multinomial")])) - 5L)
})

test_that("every roster link brms refuses for a family is refused by name", {
  sets <- frmtmb:::brms_mu_links
  roster <- names(frmtmb:::frmtmb_links)
  refused <- 0L
  for (fam in setdiff(names(sets), "multinomial")) {
    ctor <- famlink_ctor(fam)
    for (lk in setdiff(roster, sets[[fam]])) {
      expect_error(ctor(link = lk),
                   paste0("'", lk, "' is not a supported link for family '"),
                   fixed = TRUE, label = paste(fam, lk))
      refused <- refused + 1L
    }
  }
  expect_gt(refused, 300L)
  # acat's other five brms links are refused for the density, not the link
  for (lk in c("probit", "probit_approx", "cloglog", "cauchit", "softit")) {
    expect_error(acat(lk), "has not written", fixed = TRUE)
  }
})

test_that("the refusal fires on the silent wrong answer it was built for", {
  set.seed(1)
  d <- data.frame(y = rbinom(60, 1, 0.4), x = rnorm(60))
  expect_error(frm(y ~ x, d, family = bernoulli("sqrt")),
               "'sqrt' is not a supported link for family 'bernoulli'")
  # and a stats family carrying a link brms refuses is refused too
  expect_error(frm(y ~ 1, data.frame(y = rpois(20, 3)),
                   family = stats::poisson("inverse")),
               "not a supported link for family 'poisson'")
})

test_that("a custom link passes, and a malformed one names the family", {
  lk <- frmtmb:::frmtmb_links[["logit"]]
  lk$name <- "my_logit"
  expect_identical(bernoulli(lk)$link, "my_logit")
  expect_error(zero_inflated_poisson(list(name = "x")),
               "given to zero_inflated_poisson()", fixed = TRUE)
})

test_that("an unquoted link resolves as brms resolves it", {
  lk <- "probit"
  expect_identical(bernoulli(lk)$link, "probit")
  expect_identical(bernoulli(link = NULL)$link, "logit")
  expect_identical(bernoulli(link = NA)$link, "logit")
  # the default is brms's: the first link of the family's set
  expect_identical(brmsfamily("inverse.gaussian")$link, "1/mu^2")
  # a bare link name is taken as a name without being evaluated: base
  # R's identity() is a function, and a symbol bound to a DIFFERENT link
  # name is still read as the name written, as brms reads it
  probit <- "logit"
  expect_identical(bernoulli(probit)$link, "probit")
  # a value that is no link at all is refused quoting what was written
  expect_error(negbinomial(2), "'2' is not a supported link")
})

# -------------------------------------------- the family object fields

test_that("link and link_<dpar> agree with links on every registry family", {
  checked <- 0L
  for (nm in names(frmtmb:::family_registry)) {
    f <- if (nm == "multinomial") multinomial(3) else {
      frmtmb:::family_registry[[nm]]()
    }
    x <- unclass(f)
    expect_true(is.character(f$link) && length(f$link) == 1L &&
                  !is.na(f$link), label = nm)
    # read from links, never stored on the list
    expect_false(any(grepl("^link($|fun$|inv$|_)", names(x))), label = nm)
    main <- x$ord_link %||% x$links$mu
    if (!identical(x$type, "categorical") && nm != "multinomial") {
      expect_identical(f$link, main$name, label = nm)
      expect_identical(f$linkinv, main$linkinv, label = nm)
    }
    secondary <- setdiff(x$dpars, c("mu", if (nm == "multinomial" ||
                                            identical(x$type, "categorical"))
                                      x$primary_dpars))
    for (dp in secondary) {
      expect_identical(f[[c("links", dp, "name")]],
                       eval(call("$", f, paste0("link_", dp))),
                       label = paste(nm, dp))
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 30L)
})

test_that("$ does not partial-match on a family object", {
  f <- beta_binomial()
  expect_identical(f$link_phi, "log")
  expect_error(f$link_p, "partial-matched `link_phi`", fixed = TRUE)
  expect_error(f$lpd, "partial-matched `lpdf`", fixed = TRUE)
  expect_null(f$link_zi)   # absent, as on a brms family
  expect_null(f$lin)       # ambiguous, as base `$` is
  # the construction that used to hand back the links list
  expect_true(is.character(student()$link))
})

test_that("the link fields follow links through every write path", {
  lk <- frmtmb:::frmtmb_links
  base_f <- brmsfamily("gaussian")
  pos <- match("links", names(unclass(base_f)))
  paths <- list(
    dollar = function(f) { f$links$mu <- lk$log; f },
    bracket2 = function(f) { f[["links"]][["mu"]] <- lk$log; f },
    single = function(f) {
      l <- f[["links"]]; l$mu <- lk$log; f["links"] <- list(l); f
    },
    numeric = function(f) {
      l <- f[["links"]]; l$mu <- lk$log; f[[pos]] <- l; f
    },
    rebuild = function(f) {
      x <- unclass(f); x$links$mu <- lk$log
      structure(x, class = "frmtmb_family")
    },
    modify = function(f) {
      l <- f[["links"]]; l$mu <- lk$log
      utils::modifyList(f, list(links = l))
    }
  )
  for (p in names(paths)) {
    f <- paths[[p]](base_f)
    expect_identical(f$link, "log", label = p)
    expect_identical(f$linkinv, lk$log$linkinv, label = p)
  }
  # a stored element under a link field name is refused when read
  f <- base_f
  f["link"] <- "probit"
  expect_error(f$link, "stores an element named `link`", fixed = TRUE)
})

test_that("a link replaced in family_finalize() reaches family(fit)", {
  set.seed(5)
  d <- data.frame(y = rgamma(50, 3, 2), x = rnorm(50))
  fam <- frmtmb_family(
    "fin", dpars = c("mu", "shape"),
    links = list(mu = "log", shape = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dgamma(y, shape = dpars[["shape"]],
                   scale = dpars[["mu"]] / dpars[["shape"]], log = TRUE)
    },
    family_finalize = function(fam, y, aterms) {
      fam[["links"]][["shape"]] <- "softplus"
      fam
    })
  expect_identical(fam$link_shape, "log")
  fit <- frm(y ~ x, d, family = fam)
  expect_identical(family(fit)$link_shape, "softplus")
})

# ------------------------------------------------------------ mixture

test_that("mixture() refuses what brms refuses, and keeps what brms keeps", {
  expect_error(mixture(gaussian, hurdle_gamma),
               "not allowed in mixture models: hurdle_gamma", fixed = TRUE)
  expect_error(mixture(poisson, multinomial(3)),
               "not allowed in mixture models: multinomial", fixed = TRUE)
  expect_error(mixture(Beta, zero_inflated_beta),
               "not allowed in mixture models: zero_inflated_beta",
               fixed = TRUE)
  expect_error(mixture(poisson, cox), "real and integer support")
  expect_error(mixture(gaussian, poisson, nmix = 2),
               "mixture() has no argument `nmix`", fixed = TRUE)
  # brms builds both of these
  expect_s3_class(mixture(poisson, bernoulli), "frmtmb_family")
  expect_s3_class(mixture(gaussian, Beta), "frmtmb_family")
  expect_s3_class(mixture(binomial, beta_binomial), "frmtmb_family")
})

# ---------------------------------------------------- the response

test_that("an unordered factor response to an ordinal family is refused", {
  set.seed(11)
  dd <- data.frame(y = factor(sample(c("high", "low", "mid"), 60, TRUE)),
                   x = rnorm(60))
  for (fn in list(cumulative, sratio, cratio, acat)) {
    expect_error(frm(y ~ x, dd, family = fn(), dry_run = "frame"),
                 "requires either positive integers or ordered factors")
  }
  dd$yo <- factor(dd$y, levels = c("low", "mid", "high"), ordered = TRUE)
  expect_s3_class(frm(yo ~ x, dd, family = cumulative(), dry_run = "frame"),
                  "frmtmb_frame")
  # a categorical response IS an unordered factor
  expect_s3_class(frm(y ~ x, dd, family = categorical(), dry_run = "frame"),
                  "frmtmb_frame")
})

test_that("a binomial-type response without trials() is refused, as in brms", {
  set.seed(9)
  db <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40), n = 1)
  msg <- "Specifying 'trials' is required for this model."
  for (fam in list(binomial(), beta_binomial(), zero_inflated_binomial(),
                   mixture(binomial, binomial))) {
    expect_error(frm(y ~ x, db, family = fam, dry_run = "frame"), msg,
                 fixed = TRUE)
  }
  expect_error(frm(y ~ x, db, family = "binomial", dry_run = "frame"), msg,
               fixed = TRUE)
  expect_s3_class(suppressMessages(
    frm(y | trials(n) ~ x, db, family = binomial(), dry_run = "frame")),
    "frmtmb_frame")
  # a mixture of binomials is read through its components for the message
  expect_message(frm(y | trials(n) ~ x, db,
                     family = mixture(binomial, binomial),
                     dry_run = "frame"),
                 "Only 2 levels detected", fixed = TRUE)
})

test_that("a two-outcome response gets brms's bernoulli message", {
  set.seed(15)
  db <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40),
                   y2 = sample(1:2, 40, TRUE), y3 = sample(1:3, 40, TRUE),
                   f2 = factor(sample(c("a", "b"), 40, TRUE)),
                   two = 2)
  bern <- paste("Only 2 levels detected so that family 'bernoulli' might",
                "be a more efficient choice.")
  fr <- function(...) frm(..., data = db, dry_run = "frame")
  expect_message(fr(y | trials(1) ~ x, family = binomial()), bern,
                 fixed = TRUE)
  expect_message(fr(y | trials(1) ~ x, family = beta_binomial()), bern,
                 fixed = TRUE)
  expect_message(fr(y | trials(1) ~ x, family = zero_inflated_binomial()),
                 bern, fixed = TRUE)
  expect_message(fr(y2 ~ x, family = cumulative()), bern, fixed = TRUE)
  expect_message(fr(f2 ~ x, family = categorical()), bern, fixed = TRUE)
  expect_no_message(fr(y | trials(two) ~ x, family = binomial()))
  expect_no_message(fr(y3 ~ x, family = cumulative()))
  expect_no_message(fr(y ~ x, family = bernoulli()))
  expect_no_message(fr(y ~ x, family = poisson()))
})

# ------------------------------------------ brmsfamily and c(family, link)

test_that("brmsfamily() follows brms's spelling rules and refuses the rest", {
  expect_identical(brmsfamily("normal")$family, "gaussian")
  expect_identical(brmsfamily("hu_poisson")$family, "hurdle_poisson")
  expect_identical(brmsfamily("Gamma", identity)$link, "identity")
  expect_identical(brmsfamily("gamma")$family, "Gamma")
  expect_error(brmsfamily(c("poisson", "log")), "single string")
  expect_error(brmsfamily("poisson", link_sigma = "log"),
               "has no argument `link_sigma`", fixed = TRUE)
  expect_identical(brmsfamily("gaussian", "softplus")$link, "softplus")
  expect_identical(brmsfamily("com_poisson")$family, "compois")
  expect_false(exists("frm_family", envir = asNamespace("frmtmb")))
})

test_that("family = c(name, link) works in frm()", {
  expect_identical(frmtmb:::as_frmtmb_family(c("weibull", "identity"))$link,
                   "identity")
  expect_error(frmtmb:::as_frmtmb_family(c("weibull", "log", "x")),
               "not a character")
  set.seed(16)
  dw <- data.frame(y = rweibull(50, 2, 1), x = rnorm(50))
  fit <- frm(y ~ x, dw, family = c("weibull", "log"))
  ref <- frm(y ~ x, dw, family = weibull())
  expect_identical(fixef(fit), fixef(ref))
})

test_that("a brms family object keeps its secondary links in frm()", {
  skip_if_not_installed("brms")
  bfam <- brms::brmsfamily("beta_binomial", "probit", link_phi = "softplus")
  f <- frmtmb:::as_frmtmb_family(bfam)
  expect_identical(f$link, "probit")
  expect_identical(f$link_phi, "softplus")
})

# ------------------------------------- an overshoot past a link's domain

test_that("a trial step past a link's domain does not warn; a NaN start does", {
  # the band design of test-numerical-robustness.R: nlminb steps below
  # eta = 0 on the 1/mu^2 link eleven times on the way to the optimum
  set.seed(11)
  n <- 400
  x <- stats::rnorm(n)
  mig <- 1 / sqrt(pmax(0.5 + 0.2 * x, 0.05))
  d <- data.frame(y = stats::rgamma(n, shape = 20, rate = 20 / mig), x = x)
  expect_no_warning(fit <- frm(y ~ x, d, family = "inverse.gaussian"))
  expect_gt(fit$opt$nonfinite_trials, 0L)
  expect_identical(fit$opt$convergence, 0L)
  # the same fit on the unmapped objective: the same optimum, and the
  # warning the mapping removes
  raw <- function(par, fn, gr, lower, upper, control) {
    stats::nlminb(par, fn, gr, control = control, lower = lower,
                  upper = upper)
  }
  expect_warning(ref <- frm(y ~ x, d, family = "inverse.gaussian",
                            control = frmtmb_control(optimizer = raw)),
                 "NA/NaN function evaluation")
  expect_identical(ref$opt$par, fit$opt$par)
  # a density that is NaN everywhere, the start included, still warns
  broken <- frmtmb_family(
    "broken", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      log(-exp(dpars[["sigma"]])) + 0 * dpars[["mu"]]
    })
  expect_warning(try(frm(y ~ x, d, family = broken), silent = TRUE),
                 "NA/NaN function evaluation")
})

test_that("nonfinite_trials counts every optimizer run, restarts included", {
  # The reference count comes from the same nlminb path run through the
  # custom-optimizer hook, which optimize_obj() also restarts, so it sums
  # every run by construction. Mapping NaN to +Inf leaves nlminb's
  # iterates unchanged (see the test above), so the two must agree.
  counted <- 0L
  tally_opt <- function(par, fn, gr, lower, upper, control) {
    n_eval <- 0L
    f <- function(p) {
      v <- fn(p)
      n_eval <<- n_eval + 1L
      if (n_eval > 1L && is.nan(v)) {
        counted <<- counted + 1L
        return(Inf)
      }
      v
    }
    stats::nlminb(par, f, gr, control = control, lower = lower,
                  upper = upper)
  }
  mk <- function() {
    set.seed(20260916)
    data.frame(x = stats::runif(200, -2, 2),
               g = factor(rep(1:20, length.out = 200)))
  }
  # bernoulli: the first run maps its trials and the restart that
  # replaces it maps none, which is where a per-run count read 0
  d1 <- mk()
  d1$y <- stats::rbinom(200, 1, pmin(0.98, pmax(0.02, 0.5 + 0.45 * d1$x)))
  # negbinomial with a random intercept: the restart maps a trial and is
  # NOT kept, so a per-run count missed that one
  d2 <- mk()
  d2$y <- stats::rnbinom(200, size = 3, mu = pmax(
    0.1, 2 + 1.5 * d2$x + stats::rnorm(20, 0, 0.4)[d2$g]))
  cases <- list(list(y ~ x, d1, bernoulli("identity")),
                list(y ~ x + (1 | g), d2, negbinomial("identity")))
  for (cs in cases) {
    fit <- suppressWarnings(frm(cs[[1]], cs[[2]], family = cs[[3]]))
    counted <- 0L
    ref <- suppressWarnings(frm(cs[[1]], cs[[2]], family = cs[[3]],
                                control = frmtmb_control(
                                  optimizer = tally_opt)))
    expect_identical(ref$opt$par, fit$opt$par)
    expect_gt(counted, 0L)
    expect_identical(fit$opt$nonfinite_trials, counted)
  }
  # where a user looks: diagnose() and the verbose trace
  expect_output(diagnose(fit),
                paste0("Non-finite objective at ", counted, " trial points"),
                fixed = TRUE)
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(frm(y ~ x, d1, family = bernoulli("identity"),
                         verbose = TRUE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_match(grep("done", msgs, value = TRUE), "non-finite trial")
  # optim does not count them, and says nothing rather than zero
  fo <- suppressWarnings(frm(y ~ x, d1, family = poisson(),
                             control = frmtmb_control(optimizer = "optim")))
  expect_null(fo$opt$nonfinite_trials)
  expect_null(diagnose(fo, quiet = TRUE)$nonfinite_trials)
})

test_that("get_prior() answers without trials() where brms does", {
  # brms's get_prior() accepts all of these (dev/famlink-rev2-trials.R);
  # frm() refuses each, because the fit needs trials and the prior
  # table does not
  set.seed(20260916)
  n <- 60
  d <- data.frame(x = stats::rnorm(n), n = sample(4:12, n, TRUE))
  d$s <- stats::rbinom(n, d$n, 0.4)
  d$y01 <- stats::rbinom(n, 1, 0.4)
  d$one <- 1L
  d$Y <- t(sapply(d$n, function(k) stats::rmultinom(1, k, c(0.2, 0.3, 0.5))))
  cases <- list(
    list(y01 ~ x, binomial()),
    list(s ~ x, binomial()),
    list(y01 ~ x, zero_inflated_binomial()),
    list(y01 ~ 1, mixture(binomial(), binomial())),
    list(Y ~ x, multinomial(3)),
    list(Y | trials(one) ~ x, multinomial(3))
  )
  for (cs in cases) {
    expect_s3_class(get_prior(cs[[1]], data = d, family = cs[[2]]),
                    "frmtmb_prior_rows")
    expect_error(suppressMessages(frm(cs[[1]], data = d, family = cs[[2]])),
                 "trials")
  }
  # the same table as with trials written, since no slot depends on them
  expect_identical(get_prior(y01 ~ x, data = d, family = binomial()),
                   get_prior(y01 | trials(1) ~ x, data = d,
                             family = binomial()))
})
