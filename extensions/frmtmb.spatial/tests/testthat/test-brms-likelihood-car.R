# car() against brms 2.23.0, moved from frmtmb's
# tests/testthat/test-brms-likelihood.R (rows 19, 19c-esicar and 19c),
# unchanged. helper-brms.R is core's, copied.

brms_car_data <- function() {
  w <- matrix(0, 16, 16)
  g <- expand.grid(r = 1:4, c = 1:4)
  for (i in 1:16) {
    for (j in 1:16) {
      if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) w[i, j] <- 1
    }
  }
  dimnames(w) <- list(paste0("L", 1:16), paste0("L", 1:16))
  set.seed(42)
  kmat <- diag(rowSums(w)) - w + matrix(1 / (1e-3 * 16)^2, 16, 16)
  phi <- 1.2 * drop(crossprod(chol(solve(kmat)), rnorm(16)))
  loc <- factor(rep(rownames(w), each = 6), levels = rownames(w))
  d <- data.frame(loc = loc, x = rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    rnorm(nrow(d), 0, 0.5)
  list(d = d, W = w)
}

test_that("check C: row 19, car(escar) and car(icar)", {
  skip_unless_brms_fit()

  s <- brms_car_data()
  d19 <- s$d
  wmat <- s$W

  # escar, the proper CAR. brms declares the field itself, so the map is
  # the identity and the only difference is the normalizer: brms drops
  # -Nloc/2 log(2 pi) and 0.5 log det D from a density that is otherwise
  # frmtmb's, term for term.
  bform <- brms::bf(y ~ x + car(W, gr = loc, type = "escar"))
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "escar")) + gaussian(),
             data = d19, data2 = list(W = wmat))
  prior <- brms_flat_prior(bform, data = d19, family = gaussian(),
                           data2 = list(W = wmat))
  code <- brms::make_stancode(bform, data = d19, family = gaussian(),
                              prior = prior, data2 = list(W = wmat))
  sdat <- brms_standata(bform, data = d19, family = gaussian(),
                        prior = prior, data2 = list(W = wmat))
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)
  bk <- brms_car_block(fit)
  expect_lt(abs(pars[["sdcar"]] - exp(fit$estimates$theta[[1]])), 1e-14)
  expect_lt(abs(pars[["car"]] - car_rho(fit$estimates$theta[[2]])), 1e-14)
  expect_lt(max(abs(pars[["rcar"]] - brms_block_b(fit, bk))), 1e-14)
  # brms's Nneigh is the degree vector, which is what its dropped
  # normalizer is built from
  expect_identical(as.numeric(sdat$Nneigh), as.numeric(bk$aux$deg))
  brms_lp_check(bform, gaussian(), d19, fit, joint = TRUE,
                const = brms_car_const(sdat, "escar"),
                data2 = list(W = wmat))

  # icar, the intrinsic one. brms standardizes the field by sdcar and
  # adds a soft sum-to-zero term whose precision is exactly the rank-one
  # term frmtmb folds into its precision matrix, so the quadratic forms
  # agree and only the normalizer differs.
  bform2 <- brms::bf(y ~ x + car(W, gr = loc, type = "icar"))
  fit2 <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
              data = d19, data2 = list(W = wmat))
  prior2 <- brms_flat_prior(bform2, data = d19, family = gaussian(),
                            data2 = list(W = wmat))
  code2 <- brms::make_stancode(bform2, data = d19, family = gaussian(),
                               prior = prior2, data2 = list(W = wmat))
  sdat2 <- brms_standata(bform2, data = d19, family = gaussian(),
                         prior = prior2, data2 = list(W = wmat))
  sf2 <- suppressMessages(rstan::sampling(brms_stan_model(code2),
                                          data = sdat2, chains = 0))
  pars2 <- stan_pars_from_fit(fit2, sdat2, code2)
  expect_setequal(names(pars2), brms_stan_par_names(code2))
  expect_par_roundtrip(sf2, pars2)
  # brms's soft constraint and frmtmb's rank-one term are the same
  # number: normal(sum(zcar) | 0, 0.001 * Nloc) against 1 / (con_sd n)^2
  expect_match(code2, "normal_lpdf(sum(zcar) | 0, 0.001 * Nloc)",
               fixed = TRUE)
  expect_identical(brms_car_block(fit2)$aux$con_sd, 0.001)
  brms_lp_check(bform2, gaussian(), d19, fit2, joint = TRUE,
                const = brms_car_const(sdat2, "icar"),
                data2 = list(W = wmat))
})

test_that("row 19c-esicar: car(esicar) is the same model in both packages", {
  skip_unless_brms_fit()

  # esicar was an EXEMPTION here through 0.51.0, because frmtmb's
  # esicar was its icar and brms's is not. brms imposes the sum-to-zero
  # constraint exactly, and now so does frmtmb, so the row is an
  # identity like 19a and 19b.
  #
  # The two parameterizations are still not the same vector. brms
  # declares Nloc - 1 free values ON THE FIELD and derives the last as
  # minus their sum; frmtmb keeps Nloc coefficients and centers them in
  # expand_b(), which leaves the component mean an inert coordinate
  # with a tau-free density. The map is therefore a truncation of the
  # CENTERED vector and carries no Jacobian. icar is the contrast:
  # there both sides standardize by sdcar and the map does carry one.
  s <- brms_car_data()
  d19 <- s$d
  wmat <- s$W

  bform <- brms::bf(y ~ x + car(W, gr = loc, type = "esicar"))
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
             data = d19, data2 = list(W = wmat))
  prior <- brms_flat_prior(bform, data = d19, family = gaussian(),
                           data2 = list(W = wmat))
  code <- brms::make_stancode(bform, data = d19, family = gaussian(),
                              prior = prior, data2 = list(W = wmat))
  sdat <- brms_standata(bform, data = d19, family = gaussian(),
                        prior = prior, data2 = list(W = wmat))
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code),
                                         data = sdat, chains = 0))

  # brms's program, read rather than remembered: Nloc - 1 free values,
  # the last determined, the field itself in the predictor (no sdcar
  # factor, which is what separates esicar from icar), and the
  # pseudo-determinant (Nloc - 1) log tau as the normalizer.
  expect_true("zcar" %in% brms_stan_par_names(code))
  expect_match(code, "vector[Nloc - 1] zcar;", fixed = TRUE)
  expect_match(code, "rcar[Nloc] = - sum(zcar);", fixed = TRUE)
  expect_match(code, "mu[n] += rcar[Jloc[n]];", fixed = TRUE)
  expect_match(code, "(Nloc - 1) * log(tau)", fixed = TRUE)

  # frmtmb's layout did NOT change: Nloc coefficients, level-major,
  # which is what keeps ranef(), frm_linpred() and the importance layout
  # reading one value per location.
  bk <- brms_car_block(fit)
  expect_identical(bk$aux$type, "esicar")
  expect_length(bk$b_idx, 16L)
  expect_length(bk$c_idx, 16L)

  # the constraint is EXACT, not tight: the field sums to zero to
  # machine precision, where icar's residual is ~1e-9
  fld <- brms_car_field(fit, bk)
  expect_lt(abs(sum(fld)), 1e-12)

  pars <- stan_pars_from_fit(fit, sdat, code)
  expect_setequal(names(pars), brms_stan_par_names(code))
  expect_par_roundtrip(sf, pars)
  expect_length(pars[["zcar"]], 15L)
  # brms rebuilds the dropped value; it has to come back as frmtmb's
  expect_lt(abs(-sum(pars[["zcar"]]) - fld[[16]]), 1e-12)
  # no Jacobian, because neither side standardizes the field here
  expect_identical(attr(pars, "logJ"), 0)
  expect_lt(abs(pars[["sdcar"]] - exp(fit$estimates$theta[[1]])), 1e-14)

  brms_lp_check(bform, gaussian(), d19, fit, joint = TRUE,
                const = brms_car_const(sdat, "esicar"),
                data2 = list(W = wmat))
})

test_that("row 19c: bym2 keeps a different latent set", {
  skip_unless_brms()

  # THE REMAINING EXEMPTION on this row, and it is not arithmetic.
  # brms keeps the spatial and the non-spatial parts as SEPARATE latent
  # vectors, 2 * Nloc of them plus rhocar, and frmtmb integrates the
  # mixture into one dense marginal covariance over Nloc values. The
  # two joint densities are functions of different arguments, so no
  # parameter map exists and the translator refuses by name.
  s <- brms_car_data()
  d19 <- s$d
  wmat <- s$W

  bform <- brms::bf(y ~ x + car(W, gr = loc, type = "bym2"))
  prior <- brms_flat_prior(bform, data = d19, family = gaussian(),
                           data2 = list(W = wmat))
  code <- brms::make_stancode(bform, data = d19, family = gaussian(),
                              prior = prior, data2 = list(W = wmat))
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "bym2")) + gaussian(),
             data = d19, data2 = list(W = wmat))
  bk <- brms_car_block(fit)
  expect_identical(bk$aux$type, "bym2")
  expect_length(bk$b_idx, 16L)
  expect_true(all(c("zcar", "nszcar", "rhocar") %in%
                    brms_stan_par_names(code)))
  expect_error(stan_pars_from_fit(fit, brms_standata(
    bform, data = d19, family = gaussian(), prior = prior,
    data2 = list(W = wmat)), code), "different set of latent")
})

