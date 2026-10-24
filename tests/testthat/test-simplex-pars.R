# frmtmb_family(simplex_pars =): a family's extra parameter declared a
# simplex gets what a mo() simplex gets, the chart and the search off a
# face, the sampler's softmax chart, brms's draw names and a Dirichlet
# prior class of its own name. frmtmb.survival's cox() declares its
# baseline this way (c(sbhaz = "sbhaz_raw")); the cox() versions of
# these tests are in that package's test-cox-chart.R.
#
# The family here is a location mixture over three fixed centers with
# simplex weights, which is concave in the weights at a fixed location,
# as a cox() baseline is. Its centers reach the density through
# aterm_data =, the seam cox() builds its spline bases with.

tsimplex_family <- function(centers = c(-3, 0, 3)) {
  frmtmb_family(
    "tsimplex",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      w <- mo_simplex(extra$tw_raw, aterms[["tw_chart"]])
      cen <- aterms[["tw_centers"]]
      f <- 0
      for (j in seq_along(cen)) {
        f <- f + w[j] * exp(-0.5 * (y - dpars[["mu"]] - cen[j])^2)
      }
      log(f) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) 0),
    extra_pars = function(y, aterms) {
      list(tw_raw = rep(0, length(aterms[["tw_centers"]]) - 1L))
    },
    aterm_data = function(y, aterms) list(tw_centers = centers),
    simplex_pars = c(tw = "tw_raw"))
}

# weights (0.6, 0.4, 0): the third center carries none, so the maximum
# puts a weight at a face
tsimplex_data <- function(seed = 5, n = 300) {
  set.seed(seed)
  k <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  data.frame(y = 0.5 + c(-3, 0)[k] + stats::rnorm(n))
}

# The maximum by another route: EM over the weights and the shared
# location, which rises at every step
tsimplex_em_loglik <- function(y, centers = c(-3, 0, 3)) {
  w <- rep(1 / 3, 3)
  mu <- 0
  for (it in 1:5000) {
    d <- vapply(centers, function(cj) stats::dnorm(y - mu - cj), y)
    r <- sweep(d, 2L, w, "*")
    r <- r / rowSums(r)
    w <- colMeans(r)
    mu <- sum(r * outer(y, centers, "-")) / length(y)
  }
  d <- vapply(centers, function(cj) stats::dnorm(y - mu - cj), y)
  sum(log(as.vector(d %*% w)))
}

test_that("simplex_pars is validated where the family is built", {
  expect_error(tsimplex_family()$simplex_pars, NA)
  mk <- function(sp, ep = function(y, aterms) list(tw_raw = 0)) {
    frmtmb_family("t", dpars = "mu", links = list(mu = "identity"),
                  lpdf = function(y, dpars, aterms) y, extra_pars = ep,
                  simplex_pars = sp)
  }
  expect_error(mk("tw_raw"), "named")
  expect_error(mk(c(tw = NA_character_)), "named")
  expect_error(mk(c(tw = "a", tw = "b")), "named")
  # a class frmtmb already gives another meaning, and the one exception
  expect_error(mk(c(b = "tw_raw")), "prior class frmtmb already gives")
  # a dpar's name, which set_prior(class = ) reads as that dpar's
  expect_error(mk(c(mu = "tw_raw")), "prior class frmtmb already gives")
  expect_identical(mk(c(sbhaz = "tw_raw"))$simplex_pars,
                   c(sbhaz = "tw_raw"))
  expect_error(mk(c(tw = "tw_raw"), ep = NULL), "no extra_pars")
  expect_error(frmtmb_family("t", dpars = "mu",
                             links = list(mu = "identity"),
                             lpdf = function(y, dpars, aterms) y,
                             aterm_data = 1),
               "aterm_data")
  # a family that declares neither carries neither slot, so the slot
  # names every other family has are unchanged
  g <- frmtmb_family("t", dpars = "mu", links = list(mu = "identity"),
                     lpdf = function(y, dpars, aterms) y)
  expect_false(any(c("aterm_data", "simplex_pars") %in% names(unclass(g))))
})

test_that("a declared simplex at a face reaches the maximum, silently", {
  dd <- tsimplex_data()
  w <- character()
  fit <- withCallingHandlers(
    frm(y ~ 1, family = tsimplex_family(), data = dd),
    warning = function(z) {
      w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")
    })
  expect_identical(w, character())
  ref <- tsimplex_em_loglik(dd$y)
  gap <- ref - as.numeric(logLik(fit))
  expect_lt(gap, 1e-8 * abs(ref))
  # aterm_data reached the density: the third center carries nothing
  wt <- mo_simplex(fit$estimates$tw_raw)
  expect_lt(wt[3], 1e-6 * max(wt))
})

test_that("a face the likelihood rises from is left; the maximum is not", {
  # moved from test-cox-chart.R with the seam: the same guard, on a
  # family of core's own tests rather than on cox()
  dd <- tsimplex_data()
  fit <- frm(y ~ 1, family = tsimplex_family(), data = dd)
  p <- fit$opt$par
  zi <- which(names(p) == "tw_raw")
  terms <- mo_search_terms(fit$obj, fit$frame)
  expect_length(terms, 1L)
  expect_identical(terms[[1L]]$z, zi)
  expect_true(isTRUE(terms[[1L]]$face))
  tm <- terms[[1L]]
  # the guard absent: at this fit no face is one to leave
  expect_null(simplex_face_leave(fit$obj,
                                 list(par = p,
                                      objective = fit$opt$objective),
                                 tm))
  # put the largest weight at a face, where the likelihood rises from it
  w <- mo_simplex(p[zi])
  k <- which.max(w)
  w0 <- w
  w0[k] <- 0
  w0 <- w0 / sum(w0)
  p0 <- p
  p0[zi] <- mo_coords(w0)
  start <- list(par = p0, objective = fit$obj$fn(p0))
  expect_false(is.null(simplex_face_leave(fit$obj, start, tm)))
})

test_that("a declared simplex is charted and named as brms names it", {
  dd <- tsimplex_data()
  fit <- frm(y ~ 1, family = tsimplex_family(), data = dd)
  terms <- simplex_frame_terms(fit)
  expect_length(terms, 1L)
  expect_identical(terms[[1L]]$names, paste0("tw[", 1:3, "]"))
  expect_identical(terms[[1L]]$kind, "tsimplex")
  expect_identical(terms[[1L]]$D, 3L)
  fr <- mo_chart_frame(fit$frame, "softmax")
  expect_identical(fr$aterm_values[[1L]]$tw_chart, "softmax")
  expect_null(fit$frame$aterm_values[[1L]]$tw_chart)
  # the same model in other coordinates: equal objectives at the same
  # weights
  w1 <- c(0.5, 0.3, 0.2)
  obj2 <- RTMB::MakeADFun(build_objective(fr),
                          utils::modifyList(fit$estimates, list(
                            tw_raw = mo_coords(w1, "softmax"))),
                          silent = TRUE)
  p <- fit$opt$par
  p[names(p) == "tw_raw"] <- mo_coords(w1)
  v1 <- fit$obj$fn(p)
  v2 <- obj2$fn(obj2$par)
  expect_lt(abs(v1 - v2), 1e-10 * abs(v1))
})

test_that("a declared simplex takes a Dirichlet prior of its own class", {
  dd <- tsimplex_data()
  tab <- get_prior(y ~ 1, data = dd, family = tsimplex_family())
  expect_true("tw" %in% tab$class)
  fm <- frm(y ~ 1, family = tsimplex_family(), data = dd,
            prior = set_prior("dirichlet(2)", class = "tw"))
  wt <- mo_simplex(fm$estimates$tw_raw)
  expect_gt(min(wt), 1e-3 * max(wt))
  # brms's rules for a simplex, applied at the fit for a class only the
  # fit can tell is one
  expect_error(frm(y ~ 1, family = tsimplex_family(), data = dd,
                   prior = set_prior("normal(0, 1)", class = "tw")),
               "only valid prior for simplex")
  expect_error(frm(y ~ 1, family = tsimplex_family(), data = dd,
                   prior = set_prior("dirichlet(2)", class = "tw",
                                     coef = "a")),
               "nothing it can narrow to")
  expect_error(frm(y ~ 1, family = tsimplex_family(), data = dd,
                   prior = set_prior("dirichlet(c(1, 2))", class = "tw")),
               "Expected input of length 3")
  # a Dirichlet on a class that is no simplex of the model
  expect_error(frm(y ~ 1, family = gaussian(), data = dd,
                   prior = set_prior("dirichlet(2)", class = "sigma")),
               "dirichlet() is the density of a simplex", fixed = TRUE)
  expect_error(frm(y ~ 1, family = gaussian(), data = dd,
                   prior = set_prior("dirichlet(2)", class = "sbhaz")),
               "No simplex matches")
  # frmtmb's own classes are still refused where the row is written
  expect_error(set_prior("dirichlet(2)", class = "b"), "class = \"simo\"")
})

test_that("a simplex_pars value that names no extra parameter is refused", {
  # review finding 4: the declaration is checked against what
  # extra_pars() returns, which exists only once the data are in hand
  dd <- tsimplex_data()
  typo <- tsimplex_family()
  typo[["simplex_pars"]] <- c(tw = "tw_typo")
  expect_error(frm(y ~ 1, family = typo, data = dd),
               "gives no parameter called 'tw_typo'. It gives 'tw_raw'",
               fixed = TRUE)
  expect_error(get_prior(y ~ 1, data = dd, family = typo), "tw_typo",
               fixed = TRUE)
  # a simplex of one weight has no coordinate to estimate
  one <- tsimplex_family(centers = 0)
  expect_error(frm(y ~ 1, family = one, data = dd),
               "gives it no coordinate", fixed = TRUE)
  # the guard absent: the declaration that names the parameter fits
  ok <- frm(y ~ 1, family = tsimplex_family(), data = dd)
  expect_length(ok$estimates$tw_raw, 2L)
})
