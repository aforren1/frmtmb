# posterior_linpred(incl_thres = TRUE) on draws of an ordinal fit.
# brms 2.23.0 returns a draws x observations x thresholds array whose
# layer k is the argument its ordinal densities hand the link,
# disc * (thres_k - mu) for cumulative and sratio and disc * (mu -
# thres_k) for cratio and acat, with dimnames list(draw numbers, NULL,
# threshold numbers), NA past the thresholds of a thres(gr = ) level
# that has fewer (dev/fixes-thres-brms.R, at frmtmb's estimates). Up to
# frmtmb.sample 0.15.0 the argument was refused. Each draw is checked
# against brms's own R-side densities at identity link, evaluated at
# that draw's stored brms-named columns.

it_data <- function(seed = 62, n = 120) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  # level b never reaches the top category: one threshold fewer
  d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
  d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
  d
}

it_draws <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      d <- it_data()
      fit <- switch(which,
        cum = frm(y ~ x, family = cumulative(), data = d),
        sr_cs = frm(bf(y ~ x + cs(z)), family = sratio(), data = d),
        cr_disc = frm(bf(y ~ x, disc ~ 0 + z), family = cratio(), data = d),
        acat_eq = frm(y ~ x, family = acat(threshold = "equidistant"),
                      data = d),
        sr_stz = frm(y ~ x, family = sratio(threshold = "sum_to_zero"),
                     data = d),
        cum_gr = frm(yg | thres(gr = g) ~ x, family = cumulative(),
                     data = d),
        hurdle = frm(yh ~ x, family = hurdle_cumulative(), data = d))
      ds <- suppressMessages(suppressWarnings(
        frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 5)))
      cache[[which]] <<- list(d = d, fit = fit, ds = ds)
    }
    cache[[which]]
  }
})

# brms's value for draw s, from its density helper at identity link
it_brms_draw <- function(ds, s, fam, K1, thr_names, rows_K1 = NULL) {
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  f <- frmtmb.sample:::draws_fit_at(ds, s, idx)
  m <- as.matrix(ds, variable = thr_names)[s, ]
  eta <- as.numeric(frm_linpred(f, type = "link", dpar = "mu"))
  # a cs() term moves eta per threshold, as brms's mu does there
  eta_k <- frmtmb:::ord_linear_per_threshold(f, eta, NULL, NULL)
  if (!is.matrix(eta_k)) eta_k <- matrix(eta, length(eta), K1)
  disc <- as.numeric(frm_linpred(f, type = "response", dpar = "disc"))
  dens <- get(paste0("d", fam), asNamespace("brms"))
  n <- length(eta)
  out <- matrix(NA_real_, n, K1)
  for (i in seq_len(n)) {
    th <- if (is.null(rows_K1)) m else rows_K1(m, i)
    k <- length(th)
    out[i, seq_len(k)] <- dens(seq_len(k), eta = eta_k[i, seq_len(k)],
                               thres = matrix(th, 1L, k), disc = disc[i],
                               link = "identity")
  }
  out
}

it_check <- function(which, fam, K1, rows_K1 = NULL, thr_names = NULL) {
  cs <- it_draws(which)
  thr_names <- thr_names %||% paste0("b_Intercept[", seq_len(K1), "]")
  pl <- posterior_linpred(cs$ds, incl_thres = TRUE, ndraws = 4)
  expect_identical(dim(pl), c(4L, nrow(cs$d), as.integer(K1)))
  expect_identical(dimnames(pl), list(as.character(1:4), NULL,
                                      as.character(seq_len(K1))))
  rows <- frmtmb.sample:::draws_subsample(cs$ds, 4, NULL)
  for (k in c(1L, 4L)) {
    ref <- it_brms_draw(cs$ds, rows[k], fam, K1, thr_names, rows_K1)
    expect_identical(unname(is.na(pl[k, , ])), is.na(ref), label = which)
    ok <- !is.na(ref)
    # relative to the size of the predictor this draw measures
    expect_lt(max(abs(pl[k, , ][ok] - ref[ok])) / max(abs(ref[ok])),
              1e3 * .Machine$double.eps, label = which)
  }
  invisible(pl)
}

test_that("incl_thres includes the thresholds on every ordinal family", {
  it_check("cum", "cumulative", 4L)
  it_check("sr_cs", "sratio", 4L)
  it_check("cr_disc", "cratio", 4L)
  it_check("acat_eq", "acat", 4L)
  it_check("sr_stz", "sratio", 4L)
})

test_that("grouped thresholds read their level's and pad with NA", {
  cs <- it_draws("cum_gr")
  nm <- c(paste0("b_Intercept[a,", 1:4, "]"),
          paste0("b_Intercept[b,", 1:3, "]"))
  g <- cs$d$g
  pl <- it_check("cum_gr", "cumulative", 4L, thr_names = nm,
                 rows_K1 = function(m, i) {
                   if (g[i] == "a") m[1:4] else m[5:7]
                 })
  expect_true(all(is.na(pl[, g == "b", 4L])))
  expect_false(anyNA(pl[, g == "a", ]))
})

test_that("incl_thres is ignored where brms ignores it", {
  cs <- it_draws("cum")
  base <- posterior_linpred(cs$ds, ndraws = 3)
  expect_identical(posterior_linpred(cs$ds, incl_thres = TRUE, ndraws = 3,
                                     dpar = "mu"), base)
  expect_identical(posterior_linpred(cs$ds, transform = TRUE,
                                     incl_thres = TRUE, ndraws = 3),
                   posterior_linpred(cs$ds, transform = TRUE, ndraws = 3))
  expect_identical(posterior_linpred(cs$ds, incl_thres = FALSE,
                                     ndraws = 3), base)
  expect_error(posterior_linpred(cs$ds, incl_thres = NA), "incl_thres")
})

test_that("a hurdle ordinal's incl_thres is refused with the reason", {
  cs <- it_draws("hurdle")
  expect_error(posterior_linpred(cs$ds, incl_thres = TRUE, ndraws = 2),
               "linear predictor of nothing")
})
