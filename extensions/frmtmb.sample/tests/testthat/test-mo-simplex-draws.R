# A mo() simplex on the sampler (lane optima, punch round 1). The fit
# holds the simplex in frmtmb's maximum-likelihood chart, whose folds
# wall a density off into sheets; frm_sample() samples the softmax under
# brms's default dirichlet(1) instead, and reports brms's simo_ weights.
# dev/optima-findings.md has the comparison with brms 2.23.0.

mo_weak_data <- function() {
  set.seed(1)
  x <- sample(0:3, 100, TRUE)
  data.frame(x = x, y = 0.05 * x + rnorm(100))
}

mo_strong_data <- function(age = FALSE) {
  set.seed(1234)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  d <- data.frame(income = income,
                  ls = c(30, 60, 70, 75)[income] + rnorm(100, sd = 7))
  # drawn after the rest, so the other columns stay as they were
  if (age) d$age <- rnorm(100, mean = 40, sd = 10)
  d
}

test_that("a mo() simplex is drawn as brms's simo_ weights", {
  skip_on_cran()
  skip_sampler()
  fit <- frm(ls ~ mo(income), data = mo_strong_data())
  # short chains: rstan's ESS warnings are about their length
  ds <- allow_warnings(suppressMessages(frm_sample(
    fit, chains = 2, iter = 600, seed = 1, refresh = 0)),
    "Effective Samples Size")
  nm <- paste0("simo_moincome1[", 1:3, "]")
  expect_true(all(nm %in% variables(ds)))
  expect_false(any(grepl("^zeta", variables(ds))))
  W <- as.matrix(ds, variable = nm)
  # rows of the simplex, at a ratio to their own size
  expect_lt(max(abs(rowSums(W) - 1)) / max(W), 1e-12)
  # every reader hands the draws back to the model in the fit's chart:
  # the expected value of the response at a category is finite and moves
  # with the weights
  lev <- levels(mo_strong_data()$income)
  nd <- data.frame(income = factor(lev, levels = lev, ordered = TRUE))
  ep <- posterior_epred(ds, newdata = nd)
  expect_true(all(is.finite(ep)))
  expect_true(all(apply(ep, 1L, function(v) all(diff(v) > 0))))
  im <- frmtmb.sample:::draws_internal_matrix(ds)
  # the inverse maps the values back and keeps the column names
  z <- im[, nm[1:2], drop = FALSE]
  back <- t(apply(z, 1L, frmtmb:::mo_simplex))
  expect_equal(unname(back), unname(W[, 1:3]))
})

test_that("a weakly identified simplex is spread as dirichlet(1) is", {
  skip_on_cran()
  skip_sampler()
  # on 0.68.1 every draw sat at a vertex (a flat density on the softmax
  # coordinates is 1 / prod(w) on the simplex); on the maximum-likelihood
  # chart every draw sat at the barycenter with sd 0. brms 2.23.0 gives
  # each weight an sd near 0.23 here, close to dirichlet(1)'s own
  fit <- frm(y ~ mo(x), data = mo_weak_data())
  # short chains: rstan's ESS warnings are about their length
  ds <- allow_warnings(suppressMessages(frm_sample(
    fit, chains = 2, iter = 600, seed = 1, refresh = 0)),
    "Effective Samples Size")
  W <- as.matrix(ds, variable = paste0("simo_mox1[", 1:3, "]"))
  prior_sd <- sqrt((1 / 3) * (2 / 3) / 4)
  expect_true(all(apply(W, 2L, stats::sd) > 0.5 * prior_sd))
  expect_true(all(apply(W, 2L, stats::sd) < 1.5 * prior_sd))
  sp <- rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)
  expect_identical(sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0)),
                   0)
})

test_that("the default prior note names the simplex prior", {
  skip_on_cran()
  skip_sampler()
  fit <- frm(ls ~ mo(income), data = mo_strong_data())
  msg <- character()
  withCallingHandlers(
    frm_sample(fit, chains = 1, iter = 100, seed = 1, refresh = 0),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  expect_true(any(grepl("simo: dirichlet(1)", msg, fixed = TRUE)))
})

# ---- punch round 2 ----------------------------------------------------

test_that("check_laplace() compares a simplex on the weight scale", {
  skip_on_cran()
  skip_sampler()
  # the draws carry D weights and the fit D - 1 coordinates, on a sheet
  # the draws need not share: the first punch round stopped here with
  # "length(ml) == length(keep) is not TRUE"
  fit <- frm(ls ~ mo(income), data = mo_strong_data())
  msg <- character()
  cl <- withCallingHandlers(
    allow_warnings(check_laplace(fit, chains = 2, iter = 600, seed = 1,
                                 refresh = 0), "Effective Samples Size"),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  # no weight sits at 0 here, so there is no boundary to name
  expect_false(any(grepl("boundary of the simplex", msg, fixed = TRUE)))
  nm <- paste0("simo_moincome1[", 1:3, "]")
  expect_true(all(nm %in% cl$parameter))
  expect_false(any(grepl("^zeta", cl$parameter)))
  tab <- frmtmb:::summary_mo_frame(fit, 0.95)
  r <- match(nm, cl$parameter)
  expect_equal(cl$ml[r], unname(tab[, "Estimate"]))
  expect_equal(cl$wald_se[r], unname(tab[, "Est.Error"]))
  # an identified simplex: the posterior agrees with the Wald view
  expect_true(all(abs(cl$z_shift[r]) < 0.5))
  expect_true(all(cl$sd_ratio[r] > 2 / 3 & cl$sd_ratio[r] < 1.5))
})

test_that("check_laplace() names a weight at 0 instead of judging it", {
  skip_on_cran()
  skip_sampler()
  fit <- frm(y ~ mo(x), data = mo_weak_data())
  face <- attr(frmtmb:::summary_mo_frame(fit, 0.95), "face")
  expect_gt(length(face), 0)
  msg <- character()
  cl <- withCallingHandlers(
    allow_warnings(check_laplace(fit, chains = 2, iter = 600, seed = 1,
                                 refresh = 0), "Effective Samples Size"),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("boundary of the simplex", msg, fixed = TRUE)))
  expect_true(all(is.na(cl$wald_se[cl$parameter %in%
                                     paste0("simo_", face)])))
})

test_that("check_laplace() names a weight at 1 instead of judging it", {
  skip_on_cran()
  skip_sampler()
  # the interaction's simplex sits at a vertex: its weight at 1 had a
  # Wald error near 2e-4 and was flagged with sd_ratio 1262 (the
  # review's final check, c1)
  fit <- frm(ls ~ mo(income) * age, data = mo_strong_data(age = TRUE))
  vert <- "simo_moincome:age1[2]"
  msg <- character()
  cl <- withCallingHandlers(
    allow_warnings(check_laplace(fit, chains = 2, iter = 600, seed = 1,
                                 refresh = 0),
                   c("Effective Samples Size", "R-hat", "divergent",
                     "treedepth", "pairs()")),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(is.na(cl$wald_se[cl$parameter == vert]))
  q <- grep("questionable for", msg, value = TRUE, fixed = TRUE)
  expect_false(any(grepl(vert, q, fixed = TRUE)))
  b <- grep("boundary of the simplex", msg, value = TRUE, fixed = TRUE)
  expect_true(any(grepl(vert, b, fixed = TRUE)))
})

test_that("each term's simplex weights are adjacent, as brms orders them", {
  skip_on_cran()
  skip_sampler()
  set.seed(2)
  d <- data.frame(x1 = sample(0:3, 150, TRUE), x2 = sample(0:2, 150, TRUE))
  d$y <- 0.5 * d$x1 - 0.4 * d$x2 + rnorm(150)
  fit <- frm(y ~ mo(x1) + mo(x2), data = d)
  ds <- allow_warnings(suppressMessages(frm_sample(
    fit, chains = 1, iter = 300, seed = 1, refresh = 0)),
    c("Effective Samples Size", "R-hat", "divergent", "treedepth",
      "pairs()"))
  v <- colnames(ds$draws)
  s1 <- match(paste0("simo_mox11[", 1:3, "]"), v)
  s2 <- match(paste0("simo_mox21[", 1:2, "]"), v)
  expect_identical(diff(s1), rep(1L, 2))
  expect_identical(diff(s2), 1L)
  expect_identical(s2[1], s1[3] + 1L)
})

test_that("as_tmbstan() says a mo() simplex has no density there", {
  skip_on_cran()
  skip_sampler()
  fit <- frm(ls ~ mo(income), data = mo_strong_data())
  msg <- character()
  withCallingHandlers(
    suppressWarnings(as_tmbstan(fit, chains = 1, iter = 100, seed = 1,
                                refresh = 0)),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("has no density on this route", msg, fixed = TRUE)))
})
