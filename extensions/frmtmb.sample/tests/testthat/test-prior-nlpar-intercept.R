# frmtmb's spelling of a nonlinear parameter's intercept,
# set_prior(class = "Intercept", nlpar = "a"), on frm_sample(). Where
# a ~ 1 it is brms's class row class = "b", nlpar = "a", bound
# included, so the bound reaches the prior draws as well as the
# posterior. Through the first fix round it was a coefficient row whose
# bound the prior draws ignored: the prior of a went down to -1.67
# under lb = 0 (dev/reviews/2026-10-28-meme.md, S1). Where the
# parameter has more than an intercept, a bound on its intercept is a
# bound on one coefficient, which brms refuses, and so does this.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

pni_data <- function() {
  set.seed(2)
  d <- data.frame(x = stats::runif(60, 0, 3), w = stats::rnorm(60))
  d$y <- 2 * exp(-0.7 * d$x) + stats::rnorm(60, 0, 0.1)
  d
}

pni_form <- function() bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)

pni_ri <- function(fit, prior) {
  suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, prior))$ri
}

pni_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("on a ~ 1 the intercept spelling is brms's class row", {
  d <- pni_data()
  fit <- frm(pni_form(), data = d, dry_run = "objective")
  pf <- set_prior("normal(1, 1)", class = "Intercept", nlpar = "a",
                  lb = 0) +
    set_prior("normal(1, 1)", class = "Intercept", nlpar = "b", lb = 0.4)
  pb <- set_prior("normal(1, 1)", class = "b", nlpar = "a", lb = 0) +
    set_prior("normal(1, 1)", class = "b", nlpar = "b", lb = 0.4)
  rf <- pni_ri(fit, pf)
  rb <- pni_ri(fit, pb)
  # the same bounds and densities as brms's spelling (the guard-absent
  # case: class = "b", nlpar = is unchanged)
  expect_identical(rf$lower, rb$lower)
  expect_identical(rf$upper, rb$upper)
  expect_identical(lapply(rf$entries, `[[`, "dist"),
                   lapply(rb$entries, `[[`, "dist"))
  expect_identical(unname(rb$lower), c(0, 0.4))
  # the table the prior draws are read from carries the bound on the
  # class row, where brms draws prior_b_a under it
  sp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pf))
  tab <- validate_prior(sp$effective, fit)
  row <- tab[tab$class == "b" & tab$nlpar == "a" & !nzchar(tab$coef), ]
  expect_identical(as.numeric(row$lb), 0)
  expect_identical(row$prior, "normal(1, 1)")
})

test_that("the prior draws and the posterior respect the bound", {
  skip_sampler()
  d <- pni_data()
  ds <- pni_quiet(frm_sample(
    pni_form(), data = d, family = gaussian(),
    prior = set_prior("normal(1, 1)", class = "Intercept", nlpar = "a",
                      lb = 0) +
      set_prior("normal(1, 1)", class = "Intercept", nlpar = "b",
                lb = 0.4) +
      set_prior("exponential(1)", class = "sigma"),
    sample_prior = "yes", chains = 2, iter = 600, seed = 3, refresh = 0))
  pd <- prior_draws(ds)
  # brms's names for a class row's prior draws
  expect_true(all(c("b_a", "b_b") %in% names(pd)))
  expect_gte(min(pd$b_a), 0)
  expect_gte(min(pd$b_b), 0.4)
  m <- as.matrix(ds)
  expect_gte(min(m[, "b_a_Intercept"]), 0)
  expect_gte(min(m[, "b_b_Intercept"]), 0.4)
  # each draw is its normal(1, 1) truncated at the bound
  tn <- function(lo) {
    function(q) {
      (stats::pnorm(q, 1, 1) - stats::pnorm(lo, 1, 1)) /
        stats::pnorm(lo, 1, 1, lower.tail = FALSE)
    }
  }
  ks <- function(v, cdf) suppressWarnings(stats::ks.test(v, cdf))$p.value
  expect_gt(ks(pd$b_a, tn(0)), 0.01)
  expect_gt(ks(pd$b_b, tn(0.4)), 0.01)
})

test_that("with a slope, a bound on the intercept is refused", {
  d <- pni_data()
  f <- bf(y ~ a * exp(-b * x), a ~ 1 + w, b ~ 1, nl = TRUE)
  fit <- frm(f, data = d, dry_run = "objective")
  msg <- "Prior argument 'coef' may not be specified when using boundaries"
  expect_error(pni_ri(fit, set_prior("normal(1, 1)", class = "Intercept",
                                     nlpar = "a", lb = 0)),
               msg, fixed = TRUE)
  # brms's own spelling of that row is refused the same way
  expect_error(pni_ri(fit, set_prior("normal(1, 1)", class = "b",
                                     coef = "Intercept", nlpar = "a",
                                     lb = 0)),
               msg, fixed = TRUE)
  # b ~ 1 beside it is still a class row, and a's class row takes a bound
  ri <- pni_ri(fit, set_prior("normal(1, 1)", class = "Intercept",
                              nlpar = "b", lb = 0.4) +
                 set_prior("normal(0, 1)", class = "b", nlpar = "a",
                           lb = -1))
  expect_identical(unname(ri$lower[grepl("^b_", names(ri$lower))]), 0.4)
  # without a bound the intercept's own density is a coefficient row
  ri2 <- pni_ri(fit, set_prior("normal(1, 1)", class = "Intercept",
                               nlpar = "a"))
  expect_length(Filter(function(e) identical(e$dist$kind, "normal"),
                       ri2$entries), 1L)
})

# Both spellings of a ~ 1's intercept, set_prior(class = "Intercept",
# nlpar = "a") and brms's class row set_prior(class = "b", nlpar = "a"),
# name one slot. frm() resolves the pair by specificity in either order:
# the intercept spelling's density, else the class row's, and on each
# side its bound, else the class row's. frm_sample() kept the later of
# the two, so the order decided the prior: under pA + pB the posterior
# mean of b_a_Intercept was 0.470, under pB + pA 4.833
# (dev/reviews/2026-10-28-meme.md, short check of fix round 2, S3).

pni_pair <- function() {
  list(pA = set_prior("normal(5, 0.1)", class = "Intercept", nlpar = "a"),
       pB = set_prior("normal(0, 0.1)", class = "b", nlpar = "a"))
}

# the resolved density and bounds of a's one coefficient
pni_a <- function(ri, fit) {
  pt <- fit$frame$par_template
  nm <- function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  }
  hit <- Filter(function(e) identical(nm(e), "a_(Intercept)"), ri$entries)
  lo <- ri$lower[names(ri$lower) == "a_(Intercept)"]
  up <- ri$upper[names(ri$upper) == "a_(Intercept)"]
  list(dist = lapply(hit, `[[`, "dist"),
       lower = if (length(lo)) unname(lo) else NA_real_,
       upper = if (length(up)) unname(up) else NA_real_)
}

test_that("both spellings on a ~ 1 resolve by specificity in either order", {
  d <- pni_data()
  fit <- frm(pni_form(), data = d, dry_run = "objective")
  p <- pni_pair()
  bnd <- function(pr, lb) {
    s <- unclass(pr)
    s[[1L]]$lb <- lb
    structure(s, class = class(pr))
  }
  core <- function(pr) {
    pni_a(pni_quiet(frmtmb:::resolve_prior_input(fit, pr)), fit)
  }
  samp <- function(pr) pni_a(pni_ri(fit, pr), fit)
  cases <- list(
    dens = list(p$pA, p$pB),
    class_lb = list(p$pA, bnd(p$pB, 0.5)),
    icpt_lb = list(bnd(p$pA, 0.5), p$pB),
    both_lb = list(bnd(p$pA, 0.5), bnd(p$pB, 0.2)),
    icpt_bounds_only = list(set_prior("", class = "Intercept", nlpar = "a",
                                      lb = 0.5), p$pB))
  for (k in names(cases)) {
    a <- cases[[k]][[1L]]
    b <- cases[[k]][[2L]]
    ab <- samp(a + b)
    ba <- samp(b + a)
    # the order does not decide the sampled prior
    expect_identical(ab, ba, label = paste(k, "pA + pB against pB + pA"))
    # and it is the prior frm() resolves, density and bound
    fr <- core(a + b)
    expect_identical(ab$dist, fr$dist,
                     label = paste(k, "frm_sample() against frm()"))
    expect_identical(ab$lower, fr$lower, label = paste(k, "lower bound"))
  }
  # the intercept spelling's density wins, and its bound over the class
  # row's; the class row's bound is taken where it has none
  expect_identical(samp(p$pA + p$pB)$dist, samp(p$pA)$dist)
  expect_identical(samp(bnd(p$pA, 0.5) + bnd(p$pB, 0.2))$lower, 0.5)
  expect_identical(samp(p$pA + bnd(p$pB, 0.5))$lower, 0.5)
  # the guard absent: either spelling alone is unchanged
  expect_identical(samp(p$pB)$dist[[1L]]$location, 0)
  expect_identical(samp(p$pA)$dist[[1L]]$location, 5)
})

test_that("frm() gives the same fit in either order, and as the spelling alone", {
  d <- pni_data()
  p <- pni_pair()
  est <- function(pr) {
    pni_quiet(frm(pni_form(), data = d, prior = pr))$estimates$beta
  }
  alone <- est(p$pA)
  expect_identical(est(p$pA + p$pB), alone)
  expect_identical(est(p$pB + p$pA), alone)
})

test_that("a fit's own intercept prior and the call's class row sample one prior", {
  skip_sampler()
  d <- pni_data()
  p <- pni_pair()
  fit <- pni_quiet(frm(pni_form(), data = d, prior = p$pA))
  run <- function(...) {
    m <- as.matrix(pni_quiet(frm_sample(fit, ..., chains = 1, iter = 400,
                                        seed = 3, refresh = 0)))
    m[, colnames(m) != "lp__", drop = FALSE]
  }
  alone <- run()
  # the call's class row does not displace the fit's more specific row:
  # the same draws as the fit's prior alone, where the later row stood
  # and drew b_a_Intercept near 0
  expect_identical(run(prior = p$pB), alone)
  # and in one call, in either order, from a fit without a prior
  fit0 <- pni_quiet(frm(pni_form(), data = d))
  run0 <- function(pr) {
    m <- as.matrix(pni_quiet(frm_sample(fit0, prior = pr, chains = 1,
                                        iter = 400, seed = 3,
                                        refresh = 0)))
    m[, colnames(m) != "lp__", drop = FALSE]
  }
  ab <- run0(p$pA + p$pB)
  expect_identical(run0(p$pB + p$pA), ab)
  expect_identical(run0(p$pA), ab)
})

# prior_summary() of draws lists the rows the user wrote, both spellings,
# as frmtmb.sample 0.24.0 and brms do, while the sampler uses the merged
# prior. Release fix R1 printed the one merged class row in their place:
# `normal(5, 0.1)` on the class with the class row's bound, a row nobody
# wrote, and the user's `normal(10, 1)` was gone
# (dev/reviews/2026-10-29-release.md, finding 2).
test_that("prior_summary() keeps both spellings while the sampler merges them", {
  skip_sampler()
  d <- pni_data()
  pA <- set_prior("normal(5, 0.1)", class = "Intercept", nlpar = "a")
  pB2 <- set_prior("normal(10, 1)", class = "b", nlpar = "a", lb = 2)
  fit <- pni_quiet(frm(pni_form(), data = d, prior = pA))
  ds <- pni_quiet(frm_sample(fit, prior = pB2, chains = 1, iter = 100,
                             seed = 3, refresh = 0))
  tab <- as.data.frame(prior_summary(ds))
  ra <- tab[tab$nlpar == "a", c("prior", "class", "coef", "lb")]
  # the two rows as written: the coefficient row and the class row
  expect_identical(nrow(ra), 2L)
  expect_identical(sort(ra$prior), c("normal(10, 1)", "normal(5, 0.1)"))
  expect_identical(ra$coef[ra$prior == "normal(5, 0.1)"], "Intercept")
  expect_identical(ra$coef[ra$prior == "normal(10, 1)"], "")
  expect_identical(as.numeric(ra$lb[ra$prior == "normal(10, 1)"]), 2)
  # the sampler's prior is the merged one: the intercept spelling's
  # density with the class row's bound
  ri <- pni_ri(fit, pA + pB2)
  a <- pni_a(ri, fit)
  expect_identical(a$dist[[1L]]$location, 5)
  expect_identical(a$lower, 2)
  # and the printed rows, passed back, resolve to that same prior
  expect_identical(pni_a(pni_ri(fit, prior_summary(ds)), fit), a)
})
