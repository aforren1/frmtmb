# The residual correlation seam, frmtmb_register_autocor()
# (R/autocor-seam.R). A term is registered through the exported function
# only, and the registry is restored after each test that registers
# (local_frmtmb_registries()). The toy term is ?frmtmb_register_autocor's
# example, a fixed matrix like fcor(), so core's own fcor() is its
# reference: the two must agree to the bit.

myfix_term <- function() {
  list(
    build = function(spec, rows, n_all, lookup) {
      V <- lookup(spec$V, "myfix()")[rows, rows, drop = FALSE]
      L <- t(chol(V))
      list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
           logdet = sum(log(diag(L))))
    },
    whiten = function(y, mu, sigma, theta, aux) {
      u <- as.vector(aux$Linv %*% ((y - mu) / sigma))
      n <- length(y)
      lsig <- if (length(sigma) == 1L) n * log(sigma) else sum(log(sigma))
      list(u = u, logjac = -aux$logdet - lsig)
    },
    # u is a vector of innovations or a matrix of them, one column per
    # set; %*% keeps the shape
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * (aux$L %*% u)
    })
}

register_myfix <- function(...) {
  tt <- myfix_term()
  frmtmb_register_autocor("myfix", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten, ...)
}

seam_data <- function() {
  set.seed(5)
  n <- 24
  V <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + d$x + as.vector(t(chol(V)) %*% rnorm(n))
  list(d = d, V = V)
}

test_that("a registered term fits, reports and draws as core's fcor()", {
  local_frmtmb_registries()
  register_myfix()
  sd_ <- seam_data()
  d <- sd_$d
  V <- sd_$V
  f1 <- frm(y ~ x + myfix(V), data = d, data2 = list(V = V))
  f0 <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_identical(coef(f1), coef(f0))
  expect_identical(fitted(f1), fitted(f0))
  expect_identical(simulate(f1, nsim = 2, seed = 4),
                   simulate(f0, nsim = 2, seed = 4))
  expect_identical(autocor_matrix(f1), autocor_matrix(f0))
  # it joined the vocabulary as a residual correlation term, with the
  # guards every such term has
  expect_true("myfix()" %in% frm_compat_features()$name)
  expect_identical(frm_compat("myfix()", "poisson")$status, "refused")
  expect_identical(frm_compat("myfix()", "rescor")$status, "refused")
  # a predicted sigma is the term's to allow
  expect_s3_class(frm(bf(y ~ x + myfix(V), sigma ~ x), data = d,
                      data2 = list(V = V)), "frmtmb_fit")
})

test_that("sigma_varies = FALSE refuses a predicted sigma and passes one
           sigma", {
  local_frmtmb_registries()
  # the map is handed one sigma, as ?frmtmb_register_autocor says, and
  # refuses anything else
  tt <- myfix_term()
  one <- function(f) {
    function(u, mu, sigma, theta, aux) {
      if (length(sigma) != 1L) stop("sigma has ", length(sigma), " values")
      f(u, mu, sigma, theta, aux)
    }
  }
  frmtmb_register_autocor("myfix", args = "V", build = tt$build,
                          whiten = one(tt$whiten),
                          unwhiten = one(tt$unwhiten), sigma_varies = FALSE)
  sd_ <- seam_data()
  f1 <- frm(y ~ x + myfix(V), data = sd_$d, data2 = list(V = sd_$V))
  f0 <- frm(y ~ x + fcor(V), data = sd_$d, data2 = list(V = sd_$V))
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_identical(fitted(f1), fitted(f0))
  expect_identical(simulate(f1, nsim = 2, seed = 4),
                   simulate(f0, nsim = 2, seed = 4))
  expect_identical(residuals(f1, type = "pearson"),
                   residuals(f0, type = "pearson"))
  expect_error(frm(bf(y ~ x + myfix(V), sigma ~ x), data = sd_$d,
                   data2 = list(V = sd_$V)),
               "needs a constant sigma", class = "frmtmb_error")
  # the absent case
  expect_s3_class(frm(y ~ x + myfix(V), data = sd_$d,
                      data2 = list(V = sd_$V)), "frmtmb_fit")
})

test_that("a map whose two halves do not agree is refused at the frame", {
  local_frmtmb_registries()
  tt <- myfix_term()
  reg <- function(name, whiten = tt$whiten, unwhiten = tt$unwhiten) {
    frmtmb_register_autocor(name, args = "V", build = tt$build,
                            whiten = whiten, unwhiten = unwhiten)
  }
  # whiten() that gives one innovation fewer than the rows
  reg("shortu", whiten = function(y, mu, sigma, theta, aux) {
    w <- tt$whiten(y, mu, sigma, theta, aux)
    w$u <- w$u[-1L]
    w
  })
  # unwhiten() that is not the inverse of whiten(): the transpose of L
  reg("wrongun", unwhiten = function(u, mu, sigma, theta, aux) {
    mu + sigma * (t(aux$L) %*% u)
  })
  # unwhiten() that flattens a matrix of innovations
  reg("flatun", unwhiten = function(u, mu, sigma, theta, aux) {
    mu + sigma * as.vector(aux$L %*% u)
  })
  # unwhiten() that is right on one column of innovations and wrong on
  # the others
  reg("col1un", unwhiten = function(u, mu, sigma, theta, aux) {
    out <- mu + sigma * (aux$L %*% u)
    if (is.matrix(u) && ncol(u) > 1L) out[, -1L] <- 0
    out
  })
  # whiten() whose log-Jacobian is one value per row
  reg("vecjac", whiten = function(y, mu, sigma, theta, aux) {
    w <- tt$whiten(y, mu, sigma, theta, aux)
    w$logjac <- rep(w$logjac / length(y), length(y))
    w
  })
  sd_ <- seam_data()
  fit_with <- function(nm) {
    frm(stats::as.formula(paste0("y ~ x + ", nm, "(V)")), data = sd_$d,
        data2 = list(V = sd_$V))
  }
  expect_error(fit_with("shortu"),
               "'shortu' term's map is wrong: whiten() gave 23 innovation",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(fit_with("wrongun"),
               "'wrongun' term's map is wrong: whiten(unwhiten(u)) is not u",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(fit_with("flatun"),
               "unwhiten() of a matrix of innovations must give",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(fit_with("vecjac"), "must give one finite log-Jacobian",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(fit_with("col1un"),
               "unwhiten() of a matrix of innovations must give",
               fixed = TRUE, class = "frmtmb_error")
  # the absent case: core's fcor() and a correct map pass the same check
  # (frmtmb.spatial's sar() passes it in that package's test-sar.R)
  reg("rightmap")
  expect_s3_class(fit_with("rightmap"), "frmtmb_fit")
  expect_s3_class(fit_with("fcor"), "frmtmb_fit")
})

test_that("a bounded parameter takes a prior by its name", {
  local_frmtmb_registries()
  tt <- myfix_term()
  half <- function(n) seq_len(n) <= n / 2
  # the scale of the first half of the rows, on (0.5, 2)
  sc <- function(theta) 0.5 + 1.5 / (1 + exp(-theta[1]))
  frmtmb_register_autocor(
    "bscale", args = "V", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale1 = sc(theta)),
    bounds = function(aux) cbind(lower = 0.5, upper = 2),
    whiten = function(y, mu, sigma, theta, aux) {
      tt$whiten(y, mu, sigma * sc(theta)^half(length(y)), theta, aux)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      tt$unwhiten(u, mu, sigma * sc(theta)^half(NROW(u)), theta, aux)
    })
  sd_ <- seam_data()
  d <- sd_$d
  V <- sd_$V
  fit <- frm(y ~ x + bscale(V), data = d, data2 = list(V = V))
  nat1 <- function(f) {
    unname(frmtmb:::autocor_natural(f$estimates$thetaac,
                                    f$frame$autocor[[1L]]))
  }
  ap <- frmtmb_autocor_pars(fit)
  expect_identical(ap$name, "scale1")
  expect_identical(c(ap$lower, ap$upper), c(0.5, 2))
  dp <- default_prior(y ~ x + bscale(V), data = d, data2 = list(V = V))
  expect_true("scale1" %in% dp$class)
  # a tight prior on scale1 pulls it to 1.8; the same prior written on
  # the internal parameter, with the log-Jacobian of the map, is the
  # same model
  pr <- frm(y ~ x + bscale(V), data = d, data2 = list(V = V),
            prior = set_prior("normal(1.8, 0.01)", class = "scale1"))
  s1 <- nat1(pr)
  expect_lt(abs(s1 - 1.8) / 1.8, 0.01)
  # lb and ub narrow the interval
  # (the estimate, near 1, is held at the new lower end)
  pb <- allow_warnings(
    frm(y ~ x + bscale(V), data = d, data2 = list(V = V),
        prior = set_prior("", class = "scale1", lb = 1.5, ub = 1.9)),
    "a bound holds it", require = "a bound holds it")
  expect_gte(nat1(pb), 1.5)
  # an end of the interval is not inside it, as for ar
  expect_error(frm(y ~ x + bscale(V), data = d, data2 = list(V = V),
                   prior = set_prior("", class = "scale1", ub = 2)),
               "is outside the interval of bscale(V)", fixed = TRUE,
               class = "frmtmb_error")
  # a parameter whose term states no interval is refused by name, with
  # the internal spelling to use instead
  frmtmb_register_autocor(
    "uscale", args = "V", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale2 = exp(theta[1])),
    whiten = function(y, mu, sigma, theta, aux) {
      tt$whiten(y, mu, sigma * exp(theta[1] * half(length(y))), theta, aux)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(NROW(u))), theta, aux)
    })
  expect_error(frm(y ~ x + uscale(V), data = d, data2 = list(V = V),
                   prior = set_prior("normal(1, 1)", class = "scale2")),
               "whose term declares no interval", class = "frmtmb_error")
  # bounds() that do not match natural() are refused at the frame
  frmtmb_register_autocor(
    "badbound", args = "V", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale3 = exp(theta[1])),
    bounds = function(aux) cbind(lower = 0.5, upper = 2),
    whiten = function(y, mu, sigma, theta, aux) {
      tt$whiten(y, mu, sigma * exp(theta[1] * half(length(y))), theta, aux)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(NROW(u))), theta, aux)
    })
  expect_error(frm(y ~ x + badbound(V), data = d, data2 = list(V = V)),
               "natural() is not the bounded map of bounds()", fixed = TRUE,
               class = "frmtmb_error")
})

test_that("a term with a parameter reports it through natural", {
  local_frmtmb_registries()
  # sigma times a free scale on the first half of the rows: the
  # parameter is the log of that scale
  tt <- myfix_term()
  half <- function(n) seq_len(n) <= n / 2
  frmtmb_register_autocor(
    "twoscale", args = "V", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale1 = exp(theta[1])),
    whiten = function(y, mu, sigma, theta, aux) {
      n <- length(y)
      k <- exp(theta[1] * half(n))
      w <- tt$whiten(y, mu, sigma * k, theta, aux)
      w
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      n <- NROW(u)
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(n)), theta, aux)
    })
  sd_ <- seam_data()
  d <- sd_$d
  fit <- frm(y ~ x + twoscale(V), data = d, data2 = list(V = sd_$V))
  expect_true("scale1" %in% variables(fit))
  cv <- confint_varcorr(fit)
  expect_true("scale1" %in% cv$term)
  # the same model written with fcor() and a predicted sigma
  d$h <- as.numeric(half(nrow(d)))
  f2 <- frm(bf(y ~ x + fcor(V), sigma ~ h), data = d,
            data2 = list(V = sd_$V))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(f2)),
               tolerance = 1e-7)
})

test_that("a registration that would be wrong is refused and changes
           nothing", {
  local_frmtmb_registries()
  tt <- myfix_term()
  before <- frmtmb:::ac_entry_names()
  feats <- frm_compat_features()
  reg <- function(name, ...) {
    frmtmb_register_autocor(name, args = "V", build = tt$build,
                            whiten = tt$whiten, unwhiten = tt$unwhiten, ...)
  }
  ex <- function(expr, pat) {
    expect_error(expr, pat, class = "frmtmb_error")
  }
  ex(reg("fcor"), "one of frmtmb's own residual correlation terms")
  ex(reg("ar"), "one of frmtmb's own residual correlation terms")
  # sar() is frmtmb.spatial's, and only that package registers it
  ex(reg("sar", package = "otherpkg"), "only frmtmb.spatial registers it")
  ex(reg("poly"), "already means something else")
  ex(reg("gaussian"), "already means something else")
  ex(reg("us"), "already means something else")
  ex(reg("not a name"), "one syntactic name")
  ex(reg("ok1", npar = 1L), "natural =")
  ex(reg("ok1", npar = -1L), "npar")
  ex(frmtmb_register_autocor("ok1", args = "V", build = tt$build,
                             whiten = function(y) y,
                             unwhiten = tt$unwhiten), "whiten")
  ex(reg("ok1", start = NA_real_, npar = 0L), "start")
  ex(reg("ok1", bounds = 2), "bounds")
  expect_identical(frmtmb:::ac_entry_names(), before)
  expect_identical(frm_compat_features(), feats)
  # the absent case: a correct registration takes, and its own package
  # replaces it
  reg("ok1")
  expect_true("ok1" %in% frmtmb:::ac_entry_names())
  # a slot that is a function of `...` takes any arguments
  frmtmb_register_autocor("ok2", args = "V", build = tt$build,
                          whiten = function(...) tt$whiten(...),
                          unwhiten = function(...) tt$unwhiten(...))
  expect_true("ok2" %in% frmtmb:::ac_entry_names())
  reg("ok1")
  ex(reg("ok1", package = "otherpkg"), "already registered")
  reg("ok1", package = "otherpkg", replace = TRUE)
  # a covariance structure cannot take a residual term's name
  ex(frmtmb_register_covstruct("ok1", lpdf = function(b, theta, blk) 0,
                               vcov = function(theta, blk) diag(1),
                               npar = 0L),
     "cannot register 'ok1'")
})

test_that("the registrations of the tests above did not leak", {
  expect_false(any(c("myfix", "twoscale", "ok1", "ok2", "bscale",
                     "shortu") %in% frmtmb:::ac_entry_names()))
  expect_false(any(c("myfix()", "twoscale()", "ok1()") %in%
                     frm_compat_features()$name))
})

test_that("a slot that gives the wrong shape is refused at the frame", {
  local_frmtmb_registries()
  tt <- myfix_term()
  frmtmb_register_autocor("badnp", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten,
                          npar = function(aux) 2.5,
                          natural = function(theta, aux) theta)
  frmtmb_register_autocor("badst", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten,
                          npar = 1L, start = c(0, 0),
                          natural = function(theta, aux) theta)
  sd_ <- seam_data()
  expect_error(frm(y ~ x + badnp(V), data = sd_$d,
                   data2 = list(V = sd_$V)),
               "npar slot gave 2.5", class = "frmtmb_error")
  expect_error(frm(y ~ x + badst(V), data = sd_$d,
                   data2 = list(V = sd_$V)),
               "start slot gave 2 value", class = "frmtmb_error")
})

test_that("a saved fit with a registered term names the package when the
           term is gone", {
  local_frmtmb_registries()
  register_myfix(package = "frmtmbnotinstalled")
  sd_ <- seam_data()
  fit <- frm(y ~ x + myfix(V), data = sd_$d, data2 = list(V = sd_$V))
  # the registration gone, as in a session that never loaded its
  # package; local_frmtmb_registries() puts the registry back
  reg <- frmtmb:::frmtmb_autocor_registry
  reg$reg[["myfix"]] <- NULL
  expect_error(fitted(fit), "frmtmbnotinstalled",
               class = "frmtmb_extension_missing")
  expect_error(summary(fit), "frmtmbnotinstalled",
               class = "frmtmb_extension_missing")
})

test_that("frmtmb_covstruct_info() says a residual term is not a structure", {
  # the 0.75.0 release review, finding 2: "sar" was sent to the
  # structures' registry and refused as a missing package
  local_frmtmb_registries()
  register_myfix()
  for (nm in c("fcor", "sar", "myfix")) {
    e <- tryCatch(frmtmb_covstruct_info(nm), error = identity)
    expect_s3_class(e, "frmtmb_error")
    expect_false(inherits(e, "frmtmb_extension_missing"), label = nm)
    expect_match(conditionMessage(e), paste0("'", nm, "' is a residual ",
                                             "correlation term"),
                 fixed = TRUE)
    expect_match(conditionMessage(e), "frmtmb_autocor_pars(fit)",
                 fixed = TRUE)
  }
  # the guard absent: a structure and an unknown name read as before
  expect_false(frmtmb_covstruct_info("ar1")$field)
  expect_error(frmtmb_covstruct_info("nosuchname"),
               "no covariance structure 'nosuchname'", fixed = TRUE)
})
