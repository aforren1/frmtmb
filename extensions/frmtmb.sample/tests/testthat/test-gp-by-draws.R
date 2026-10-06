# gp(x, by = f) on draws: brms's draw names for every sub-GP's sd and
# length scale, on brms's scales, and a prediction at a position the fit
# did not see that is a DRAW of the field there given the draw's own
# values at the fitted positions, as brms's posterior_epred() makes it
# (brms:::.predictor_gp_new()), and not that conditional's mean.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

gpd_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(5)
      n <- 60
      d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                      f = factor(rep(c("a", "b"), length.out = n)))
      d$y <- 0.5 + ifelse(d$f == "a", sin(d$x), cos(d$x)) +
        stats::rnorm(n, 0, 0.3)
      # the formula route: from a fit of an exact gp() the chain does not
      # move on some seeds, under the same priors (dev/test-backlog.md;
      # dev/gpby-findings.md, Punch round 1, m1)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(bf(y ~ gp(x, by = f)), family = gaussian(), data = d,
                   chains = 1, iter = 300, refresh = 0, seed = 4)))
      cache <<- list(d = d, fit = ds$fit, ds = ds)
    }
    cache
  }
})

test_that("draws carry brms's sdgp_ and lscale_ on brms's scales", {
  cs <- gpd_case()
  v <- variables(cs$ds)
  expect_true(all(c("sdgp_gpxfa", "sdgp_gpxfb", "lscale_gpxfa",
                    "lscale_gpxfb") %in% v))
  expect_false(any(grepl("^theta_", v)))
  m <- as.matrix(cs$ds)
  im <- draws_internal_matrix(cs$ds)
  for (j in 1:2) {
    lev <- c("a", "b")[j]
    bk <- Filter(function(b) b$covstruct == "gp",
                 cs$fit$frame$re_blocks)[[j]]
    th <- im[, draws_index(cs$ds)$theta[bk$theta_idx], drop = FALSE]
    xr <- cs$d$x[cs$d$f == lev]
    expect_lt(max(abs(m[, paste0("sdgp_gpxf", lev)] / exp(th[, 1]) - 1)),
              1e-12)
    expect_lt(max(abs(m[, paste0("lscale_gpxf", lev)] * diff(range(xr)) /
                        exp(th[, 2]) - 1)), 1e-12)
  }
})

test_that("a prediction past the positions draws the field there", {
  cs <- gpd_case()
  # row 1 at an observed position of level a, rows 2 and 3 past the data
  nd <- data.frame(x = c(cs$d$x[1], 7, 7.5),
                   f = factor(c("a", "a", "b"), levels = c("a", "b")))
  ids <- rep(1L, 1500L)
  set.seed(1)
  e1 <- posterior_epred(cs$ds, newdata = nd, draw_ids = ids)
  # at an observed position the field is the draw's own value: nothing
  # is drawn there
  expect_identical(length(unique(e1[, 1])), 1L)
  # past the data it is drawn from the conditional law at the draw's
  # parameters, whose variance frm_lp_basis() reports
  sh <- draws_fit_at(cs$ds, 1L)
  # the same draw with the kriging draw off: its conditional law
  sh[["krige_draw"]] <- NULL
  ev <- frmtmb::frm_lp_basis(sh, newdata = nd, re_formula = NA)$extra_var
  expect_identical(ev[1], 0)
  r <- apply(e1[, 2:3], 2, stats::var) / ev[2:3]
  # four Monte Carlo standard errors of a variance from n draws, whose
  # relative standard error is sqrt(2 / (n - 1))
  nn <- nrow(e1)
  expect_lt(max(abs(r - 1)), 4 * sqrt(2 / (nn - 1)))
  # the two rows past the data are draws of different fields (two
  # levels), so they are uncorrelated: four standard errors of a
  # correlation of 0, 1 / sqrt(n)
  expect_lt(abs(stats::cor(e1[, 2], e1[, 3])), 4 / sqrt(nn))
  # and the conditional MEAN is what the estimate-at-a-draw gives
  mu <- as.numeric(frm_linpred(sh, newdata = nd, re_formula = NA,
                               type = "response"))
  expect_lt(abs(mean(e1[, 2]) - mu[2]), 4 * sqrt(ev[2] / 1500))
})

test_that("a plain gp() draws the field past its positions too", {
  # the same statement without a by variable, which 0.67.0 could fit:
  # its posterior_epred() past the data was the conditional mean at
  # every repetition of one draw, with no spread at all
  skip_sampler()
  set.seed(9)
  d <- data.frame(x = sort(stats::runif(40, 0, 5)))
  d$y <- sin(d$x) + stats::rnorm(40, 0, 0.3)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ gp(x)), family = gaussian(), data = d, chains = 1,
               iter = 300, refresh = 0, seed = 4)))
  nd <- data.frame(x = c(6, 7))
  set.seed(3)
  # a draw at the median sd
  sdg <- draws_internal_matrix(ds)[, draws_index(ds)$theta[1]]
  id <- which.min(abs(sdg - stats::median(sdg)))
  e <- posterior_epred(ds, newdata = nd, draw_ids = rep(id, 1500L))
  sh <- draws_fit_at(ds, id)
  sh[["krige_draw"]] <- NULL
  ev <- frmtmb::frm_lp_basis(sh, newdata = nd, re_formula = NA)$extra_var
  expect_true(all(ev > 0))
  nn <- nrow(e)
  expect_lt(max(abs(apply(e, 2, stats::var) / ev - 1)),
            4 * sqrt(2 / (nn - 1)))
  # two positions of ONE field covary as the conditional law says
  S <- frmtmb::frm_lp_basis(sh, newdata = nd, re_formula = NA,
                            extra_cov = TRUE)$extra_cov
  rho <- stats::cov2cor(as.matrix(S))[1, 2]
  # four standard errors of a correlation near rho, (1 - rho^2) / sqrt(n)
  expect_lt(abs(stats::cor(e[, 1], e[, 2]) - rho),
            4 * (1 - rho^2) / sqrt(nn))
})

test_that("the sampling defaults are brms's sdgp and lscale rows", {
  # brms's prior_gp(): its def_scale_prior on class sdgp, and on class
  # lscale an inv_gamma tuned to each sub-GP's distances
  # (def_lscale_prior()), which frmtmb.sample solves for itself
  skip_if_not_installed("brms")
  set.seed(1)
  n <- 40
  dd <- data.frame(x = stats::runif(n, 0, 5), z = stats::runif(n, 0, 3),
                   f = factor(sample(c("a", "b", "c"), n, TRUE)))
  dd$y <- sin(dd$x) + stats::rnorm(n, 0, 0.3)
  flat <- function(p) ifelse(is.na(p) | !nzchar(p), "(flat)", p)
  for (fm in c("y ~ gp(x, by = f)", "y ~ gp(x, z, by = f, iso = FALSE)",
               "y ~ gp(x, z, by = f, k = 5)")) {
    b <- as.data.frame(brms::default_prior(brms::bf(stats::as.formula(fm)),
                                           data = dd))
    g <- as.data.frame(get_prior(frmtmb::bf(stats::as.formula(fm)),
                                 data = dd, route = "sample"))
    keep <- c("sdgp", "lscale")
    b <- b[b$class %in% keep, c("prior", "class", "coef")]
    g <- g[g$class %in% keep, c("prior", "class", "coef")]
    m <- merge(b, g, by = c("class", "coef"), all = TRUE)
    expect_identical(nrow(m), nrow(b))
    expect_identical(flat(m$prior.x), flat(m$prior.y), info = fm)
  }
})
