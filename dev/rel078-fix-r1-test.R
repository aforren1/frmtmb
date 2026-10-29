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
