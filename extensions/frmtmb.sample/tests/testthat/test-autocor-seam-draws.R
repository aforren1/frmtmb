# The draws of a model with a residual correlation term of core's
# registry carry the name the term gives its parameter, read off the
# registry (frmtmb_autocor_pars()), and a bounded parameter gets brms's
# flat default prior on its interval. frmtmb.spatial's sar() reports
# lagsar or errorsar; the terms registered here, through the exported
# frmtmb_register_autocor() alone, report lagsar or a name frmtmb.sample
# has never seen, so the test needs no other extension.

skip_on_cran()

reg_lagterm <- function(par = "lagsar") {
  frmtmb::frmtmb_register_autocor(
    "lagtest", args = "W",
    build = function(spec, rows, n_all, lookup) {
      W <- lookup(spec$W, "lagtest()")[rows, rows, drop = FALSE]
      ev <- Re(eigen(W, only.values = TRUE)$values)
      list(W = W, ev = ev, lo = 1 / min(ev), hi = 1 / max(ev))
    },
    npar = 1L,
    start = function(aux) stats::qlogis(-aux$lo / (aux$hi - aux$lo)),
    natural = function(theta, aux) {
      stats::setNames(aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1])),
                      par)
    },
    bounds = function(aux) cbind(lower = aux$lo, upper = aux$hi),
    whiten = function(y, mu, sigma, theta, aux) {
      rho <- aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1]))
      r <- y - rho * as.vector(aux$W %*% y) - mu
      # sigma_varies = FALSE: sigma is one value
      list(u = r / sigma,
           logjac = sum(log1p(-rho * aux$ev)) - length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      rho <- aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1]))
      solve(diag(nrow(aux$W)) - rho * aux$W, mu + sigma * u)
    },
    sigma_varies = FALSE)
}

lag_case <- function() {
  set.seed(5)
  n <- 25
  P <- matrix(0, n, n)
  P[cbind(1:(n - 1), 2:n)] <- 1
  P <- P + t(P)
  W <- P / rowSums(P)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- as.vector(solve(diag(n) - 0.5 * W,
                         1 + 0.5 * d$x + stats::rnorm(n)))
  list(d = d, W = W)
}

test_that("a residual term's parameter is a derived column under its brms
           name", {
  frmtmb::frmtmb_local_registries()
  reg_lagterm()
  cs <- lag_case()
  fit <- frm(y ~ x + lagtest(W), data = cs$d, data2 = list(W = cs$W))
  dv <- frmtmb.sample:::draws_cov_cols(fit)
  expect_identical(dv$names, "lagsar")
  expect_true("thetaac" %in% names(dv$cols))
  # the absent case: a model without the term has no such column
  f0 <- frm(y ~ x, data = cs$d)
  expect_false("lagsar" %in% frmtmb.sample:::draws_cov_cols(f0)$names)
})

test_that("a parameter that is not called lagsar is a column under its own
           name", {
  frmtmb::frmtmb_local_registries()
  reg_lagterm("toyphi")
  cs <- lag_case()
  fit <- frm(y ~ x + lagtest(W), data = cs$d, data2 = list(W = cs$W))
  dv <- frmtmb.sample:::draws_cov_cols(fit)
  expect_identical(dv$names, "toyphi")
  expect_true("thetaac" %in% names(dv$cols))
  # and its default prior is the flat one on its interval, by its name
  dp <- frmtmb.sample:::default_priors_for(fit)
  # a natural-scale class is held as the dpar of an Intercept
  cls <- vapply(unclass(dp), function(s) {
    if (nzchar(s$dpar %||% "")) s$dpar else s$class
  }, "")
  expect_true("toyphi" %in% cls)
  expect_false("lagsar" %in% cls)
})

test_that("the sampled draws carry lagsar, read off each draw's internal
           parameter", {
  skip_sampler()
  frmtmb::frmtmb_local_registries()
  reg_lagterm()
  cs <- lag_case()
  fit <- frm(y ~ x + lagtest(W), data = cs$d, data2 = list(W = cs$W))
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, warmup = 150, refresh = 0,
               seed = 5)))
  m <- as.matrix(ds)
  expect_true(all(c("thetaac_1", "lagsar") %in% colnames(m)))
  aux <- fit$frame[["autocor"]][[1L]][["aux"]]
  want <- aux$lo + (aux$hi - aux$lo) / (1 + exp(-m[, "thetaac_1"]))
  expect_equal(unname(m[, "lagsar"]), unname(want), tolerance = 1e-12)
})

test_that("the default prior on a bounded parameter is flat on its
           interval, as brms's on rho", {
  skip_sampler()
  frmtmb::frmtmb_local_registries()
  # a term whose likelihood does not depend on its parameter: the
  # draws of the parameter are then its prior. brms's prior on lagsar
  # and errorsar is flat on the interval; a prior flat on the internal
  # parameter instead puts the draws at the ends.
  lo <- -0.4
  hi <- 0.9
  frmtmb::frmtmb_register_autocor(
    "flatphi", args = "V",
    build = function(spec, rows, n_all, lookup) list(n = length(rows)),
    npar = 1L, start = function(aux) 0,
    natural = function(theta, aux) {
      c(toyphi = lo + (hi - lo) / (1 + exp(-theta[1])))
    },
    bounds = function(aux) cbind(lower = lo, upper = hi),
    whiten = function(y, mu, sigma, theta, aux) {
      list(u = (y - mu) / sigma, logjac = -length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) mu + sigma * u,
    sigma_varies = FALSE)
  set.seed(3)
  d <- data.frame(x = stats::rnorm(30))
  d$y <- 1 + d$x + stats::rnorm(30)
  fit <- allow_warnings(
    frm(y ~ x + flatphi(V), data = d, data2 = list(V = 1)),
    "Standard errors are not available")
  msgs <- character()
  ds <- withCallingHandlers(
    frm_sample(fit, chains = 4, iter = 3000, warmup = 500, refresh = 0,
               seed = 11),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  # the announcement names the parameter and its prior
  ann <- grep("toyphi ", unlist(strsplit(msgs, "\n", fixed = TRUE)),
              value = TRUE, fixed = TRUE)
  expect_length(ann, 1L)
  expect_match(ann, "(flat), lb = -0.4, ub = 0.9  [flat on its interval",
               fixed = TRUE)
  z <- as.matrix(ds)[, "toyphi"]
  expect_true(all(z > lo & z < hi))
  # each decile of the interval holds a tenth of the draws, within four
  # standard errors at the draws' effective size
  ess <- posterior::ess_basic(z)
  se <- sqrt(0.1 * 0.9 / ess)
  frac <- tabulate(findInterval(z, lo + (hi - lo) * (1:9) / 10) + 1L, 10L) /
    length(z)
  expect_lt(max(abs(frac - 0.1)) / se, 4)
  # and the mean is the interval's middle, within four standard errors
  expect_lt(abs(mean(z) - (lo + hi) / 2) / (sd(z) / sqrt(ess)), 4)
})

test_that("the announced interval of a default prior is rounded, and the
           prior itself is not", {
  frmtmb::frmtmb_local_registries()
  reg_lagterm()
  cs <- lag_case()
  fit <- frm(y ~ x + lagtest(W), data = cs$d, data2 = list(W = cs$W))
  ap <- frmtmb::frmtmb_autocor_pars(fit)
  dp <- frmtmb.sample:::default_priors_for(fit)
  s <- Filter(function(e) identical(e$dpar, "lagsar"), unclass(dp))[[1L]]
  # the prior holds the interval exactly
  expect_identical(c(s$dist$lower, s$dist$upper), c(ap$lower, ap$upper))
  msgs <- character()
  withCallingHandlers(
    frmtmb.sample:::announce_default_priors(dp, NULL),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  ln <- grep("lagsar ", unlist(strsplit(msgs, "\n", fixed = TRUE)),
             value = TRUE, fixed = TRUE)
  expect_length(ln, 1L)
  # brms's (flat), with the interval as lb and ub
  shown <- regmatches(ln, regexec("[(]flat[)], lb = ([^,]+), ub = ([^ ]+)",
                                  ln))[[1L]]
  expect_length(shown, 3L)
  # six digits at most, and each within a relative 1e-5 of the end
  expect_true(all(nchar(gsub("[^0-9]", "", shown[2:3])) <= 6L))
  ends <- as.numeric(shown[2:3])
  expect_lt(max(abs(ends - c(ap$lower, ap$upper)) /
                  abs(c(ap$lower, ap$upper))), 1e-5)
})
