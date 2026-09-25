# A latent-class or hidden Markov model's location is several dpars
# (theta1, theta2, ... for lca(); mu1, mu2, ... for hmm()), and frmtmb
# addresses each by dpar, as brms does for the several-location
# families it has (categorical, mixture). Through frmtmb 0.62.0 a class
# "b" or "Intercept" prior with no dpar reached every class or state at
# once, and default_prior() listed one row with an empty dpar
# (dev/mvprior-log/extfam-base.txt, -lane.txt in the frmtmb tree, script
# dev/mvprior-extfam.R).

pd_msg <- function(expr) {
  e <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (inherits(e, "frmtmb_error")) conditionMessage(e) else ""
}

pd_reached <- function(formula, family, data, pl) {
  des <- frmtmb:::prior_design(formula, data, family, list())
  r <- frmtmb:::resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  sort(unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  })))
}

pd_lca_data <- function() {
  set.seed(3)
  n <- 200
  cl <- stats::rbinom(n, 1, 0.4) + 1
  pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
  Y <- matrix(0L, n, 4)
  for (j in 1:4) Y[, j] <- 1L + stats::rbinom(n, 1, pr[cl, j])
  d <- data.frame(x = stats::rnorm(n))
  d$Y <- Y
  d
}

pd_hmm_data <- function() {
  set.seed(11)
  do.call(rbind, lapply(1:15, function(id) {
    s <- integer(20)
    s[1] <- 1L
    for (t in 2:20) {
      s[t] <- sample.int(2, 1, prob = if (s[t - 1] == 1) c(.9, .1) else
        c(.2, .8))
    }
    data.frame(id = id, t = 1:20, x = stats::rnorm(20),
               y = stats::rnorm(20, c(0, 3)[s]))
  }))
}

test_that("lca(): b and Intercept with no dpar are refused", {
  d <- pd_lca_data()
  for (p in list(set_prior("normal(0, 1)", class = "b"),
                 set_prior("normal(0, 1)", class = "Intercept"))) {
    expect_match(pd_msg(validate_prior(p, bf(Y ~ x), data = d,
                                       family = lca(K = 3))),
                 "several distributional parameters, theta1, theta2")
  }
})

test_that("hmm(): b and Intercept with no dpar are refused", {
  d <- pd_hmm_data()
  fam <- hmm(K = 2, gaussian(), time = t, group = id)
  for (p in list(set_prior("normal(0, 1)", class = "b"),
                 set_prior("normal(0, 1)", class = "Intercept"))) {
    expect_match(pd_msg(validate_prior(p, bf(y ~ x), data = d,
                                       family = fam)),
                 "several distributional parameters, mu1, mu2")
  }
})

test_that("default_prior() lists each class or state by dpar", {
  tl <- as.data.frame(default_prior(bf(Y ~ x), data = pd_lca_data(),
                                    family = lca(K = 3)))
  expect_true(all(c("theta1", "theta2") %in% tl$dpar[tl$class == "b"]))
  expect_false(any(tl$class %in% c("b", "Intercept") & !nzchar(tl$dpar)))
  th <- as.data.frame(default_prior(bf(y ~ x), data = pd_hmm_data(),
                                    family = hmm(K = 2, gaussian(),
                                                 time = t, group = id)))
  expect_true(all(c("mu1", "mu2") %in% th$dpar[th$class == "b"]))
  expect_false(any(th$class %in% c("b", "Intercept") & !nzchar(th$dpar)))
})

# These pass on 0.62.0 too: the dpar spelling reaches one class or
# state, which is what ?lca shows.
test_that("a prior with dpar reaches that class or state only", {
  expect_identical(pd_reached(bf(Y ~ x), lca(K = 3), pd_lca_data(),
                              set_prior("normal(0, 1)", class = "b",
                                        dpar = "theta2")), "theta2_x")
  expect_identical(pd_reached(bf(y ~ x), hmm(K = 2, gaussian(), time = t,
                                             group = id), pd_hmm_data(),
                              set_prior("normal(0, 1)", class = "b",
                                        dpar = "mu1")), "mu1_x")
})
