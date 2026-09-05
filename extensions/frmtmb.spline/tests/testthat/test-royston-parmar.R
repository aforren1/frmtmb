## The Royston-Parmar family against flexsurv.
##
## The headline is an IDENTITY, not an agreement: flexsurv's fitted
## coefficients are pushed through frmtmb's objective and the two log
## likelihoods are the same number. Two packages that agree to 1e-15 on
## an arbitrary parameter point are evaluating one function.

sp_bc <- function() {
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  bc
}

#' flexsurv's parameter vector, reordered into frmtmb's.
#'
#' frmtmb puts every mu coefficient in `beta` and the other dpars'
#' intercepts in `betad`, so gamma0 leads, the covariate effects on
#' gamma0 follow, and gamma1 upward come last.
#' @noRd
sp_fs_par <- function(fs) {
  cf <- fs$res[, "est"]
  gam <- cf[grep("^gamma", names(cf))]
  c(gam[1L], cf[!grepl("^gamma", names(cf))], gam[-1L])
}

test_that("the objective at flexsurv's coefficients IS flexsurv's log likelihood", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  skip_on_cran()
  bc <- sp_bc()
  for (scale in c("hazard", "odds", "normal")) {
    for (k in c(0L, 1L, 3L)) {
      fs <- flexsurv::flexsurvspline(
        survival::Surv(recyrs, censrec) ~ group, data = bc, k = k,
        scale = scale)
      kn <- fs$knots
      fam <- royston_parmar(
        knots = if (k > 0L) kn[2:(length(kn) - 1L)] else numeric(0),
        bknots = kn[c(1L, length(kn))], scale = scale)
      o <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                       family = fam, data = bc, dry_run = "objective")
      par <- stats::setNames(sp_fs_par(fs), names(o$obj$par))
      ll <- -o$obj$fn(par)
      expect_equal(ll, fs$loglik, tolerance = 1e-10,
                   info = paste(scale, k))
      # relative, so the assertion means the same at every scale
      expect_lt(abs(ll - fs$loglik) / abs(fs$loglik), 1e-13)
      # THE REACH OF THIS TEST, asserted so that the next reader knows
      # what it does not cover. The scored log S loses accuracy when
      # -log S passes about 19.2 and is floored past 36; every censored
      # row here sits an order of magnitude below that, so this identity
      # says nothing about the region rp_floored() refuses.
      fit_reach <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                               family = fam, data = bc)
      reach <- rp_floored(fit_reach, action = "report")
      expect_lt(reach$max_nlogS, 5)
      expect_equal(reach$n_censored_floored, 0L)
    }
  }
})

test_that("the two optima agree, and frmtmb's is never the worse", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  skip_on_cran()
  bc <- sp_bc()
  for (k in c(0L, 1L)) {
    fs <- flexsurv::flexsurvspline(
      survival::Surv(recyrs, censrec) ~ group, data = bc, k = k,
      scale = "hazard")
    kn <- fs$knots
    fit <- frmtmb::frm(
      frmtmb::bf(recyrs | cens(censored) ~ group),
      family = royston_parmar(
        knots = if (k > 0L) kn[2:(length(kn) - 1L)] else numeric(0),
        bknots = kn[c(1L, length(kn))]),
      data = bc)
    expect_gte(as.numeric(stats::logLik(fit)) - fs$loglik, -1e-6)
    expect_equal(as.numeric(stats::logLik(fit)), fs$loglik, tolerance = 1e-6)
    expect_equal(unname(unlist(frmtmb::fixef(fit))),
                 unname(sp_fs_par(fs)), tolerance = 1e-2)
  }
})

test_that("the default knots are flexsurv's default knots", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  bc <- sp_bc()
  fs <- flexsurv::flexsurvspline(
    survival::Surv(recyrs, censrec) ~ group, data = bc, k = 2)
  fit <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                     family = royston_parmar(df = 3), data = bc)
  # the knots live on the finalized family the fit carries
  kn <- environment(stats::family(fit)[["lpdf"]])$allknots
  expect_equal(unname(kn), unname(fs$knots), tolerance = 1e-10)
  expect_equal(as.numeric(stats::logLik(fit)), fs$loglik, tolerance = 1e-5)
})

test_that("a covariate on gamma1 is a time-varying effect", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  skip_on_cran()
  bc <- sp_bc()
  fs <- flexsurv::flexsurvspline(
    survival::Surv(recyrs, censrec) ~ group, data = bc, k = 1,
    scale = "hazard", anc = list(gamma1 = ~ group))
  kn <- fs$knots
  fit <- frmtmb::frm(
    frmtmb::bf(recyrs | cens(censored) ~ group, gamma1 ~ group),
    family = royston_parmar(knots = kn[2:(length(kn) - 1L)],
                            bknots = kn[c(1L, length(kn))]),
    data = bc)
  # the same model, so the same maximized log likelihood
  expect_equal(as.numeric(stats::logLik(fit)), fs$loglik, tolerance = 1e-4)
  # and it is a bigger model than the proportional-hazards one
  ph <- frmtmb::frm(
    frmtmb::bf(recyrs | cens(censored) ~ group),
    family = royston_parmar(knots = kn[2:(length(kn) - 1L)],
                            bknots = kn[c(1L, length(kn))]),
    data = bc)
  expect_gt(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(ph)))
})

test_that("weights and truncation behave as the compat table claims", {
  skip_if_not_installed("flexsurv")
  bc <- sp_bc()
  base <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                      family = royston_parmar(df = 2), data = bc)
  bc$w <- 1
  fw <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) + weights(w) ~ group),
                    family = royston_parmar(df = 2), data = bc)
  expect_equal(as.numeric(stats::logLik(fw)),
               as.numeric(stats::logLik(base)), tolerance = 1e-6)
  bc$w2 <- 2
  fw2 <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) + weights(w2) ~ group),
                     family = royston_parmar(df = 2), data = bc)
  expect_equal(as.numeric(stats::logLik(fw2)),
               2 * as.numeric(stats::logLik(base)), tolerance = 1e-4)
  # truncation at a bound below every observed time barely moves the
  # fit, which is the case where an lcdf that is right must barely move
  # it: the window it divides by is almost the whole line
  expect_gt(min(bc$recyrs), 0.02)
  ftr <- frmtmb::frm(
    frmtmb::bf(recyrs | cens(censored) + trunc(lb = 0.01) ~ group),
    family = royston_parmar(df = 2), data = bc)
  expect_s3_class(ftr, "frmtmb_fit")
  expect_true(is.finite(as.numeric(stats::logLik(ftr))))
  expect_equal(unlist(frmtmb::fixef(ftr))[["mu.groupPoor"]],
               unlist(frmtmb::fixef(base))[["mu.groupPoor"]],
               tolerance = 0.05)
})

test_that("simulate draws from the density it scores", {
  skip_if_not_installed("flexsurv")
  bc <- sp_bc()
  fit <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                     family = royston_parmar(df = 2), data = bc)
  set.seed(4)
  s <- stats::simulate(fit, nsim = 5)
  expect_equal(dim(as.matrix(s)), c(nrow(bc), 5L))
  expect_true(all(as.matrix(s) > 0))
  expect_true(all(is.finite(as.matrix(s))))
  # refitting the draws recovers the coefficients they were drawn from
  d2 <- bc
  d2$recyrs <- as.matrix(s)[, 1L]
  d2$censored <- 0
  kn <- environment(stats::family(fit)[["lpdf"]])$allknots
  f2 <- frmtmb::frm(frmtmb::bf(recyrs ~ group),
                    family = royston_parmar(knots = kn[2], bknots = kn[c(1, 3)]),
                    data = d2)
  expect_equal(unlist(frmtmb::fixef(f2))[["mu.groupPoor"]],
               unlist(frmtmb::fixef(fit))[["mu.groupPoor"]],
               tolerance = 0.35)
})

test_that("the family refuses what it cannot do, by name", {
  expect_error(royston_parmar(df = 0), "at least 1")
  expect_error(royston_parmar(df = 2.5), "at least 1")
  expect_error(royston_parmar(knots = "a"), "numeric vector on the LOG")
  expect_error(royston_parmar(df = 4, knots = c(0, 1)), "disagree")
  expect_error(royston_parmar(bknots = 1), "two boundary knots")
  expect_error(royston_parmar(scale = "hazards"), "should be one of")
  # a density with no knots yet cannot be evaluated
  fam <- royston_parmar(df = 2)
  expect_error(fam[["lpdf"]](1, list(mu = 0, gamma1 = 1, gamma2 = 0), list()),
               "no knots yet")
})

test_that("a non-positive survival time is refused", {
  skip_if_not_installed("flexsurv")
  bc <- sp_bc()
  bc$recyrs[1L] <- 0
  expect_error(frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                           family = royston_parmar(df = 2), data = bc),
               "strictly positive")
})

test_that("the monotonicity floor keeps the log density finite", {
  # a spline whose derivative in log time is negative is not a hazard,
  # and the density must still be a number the optimizer can use
  fam <- royston_parmar(knots = numeric(0), bknots = c(-1, 1))
  ll <- fam[["lpdf"]](c(0.5, 1, 2), list(mu = 0, gamma1 = -1), list())
  expect_true(all(is.finite(ll)))
  expect_true(all(ll < -20))
  # and it is inert where the model is sane
  ok <- fam[["lpdf"]](c(0.5, 1, 2), list(mu = 0, gamma1 = 1), list())
  x <- log(c(0.5, 1, 2))
  expect_equal(ok, x + log(1) - x - exp(x), tolerance = 1e-14)
})
