# ---- the 0.77.0 consolidation: (flat) meets sigmainn and mointerop ------

test_that("a bounds-only (flat) sigma under ar(cov = TRUE) is on sigma itself", {
  # sigma is brms's innovation sd there, so set_prior("", class =
  # "sigma", lb = 0.1) is flat on that sd above 0.1: one flat entry,
  # reading no ar, and the bound log(0.1) on log sigma
  d <- fp_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  r <- fp_rp(fit, set_prior("", class = "sigma", lb = 0.1))
  sg <- Filter(function(e) identical(e$comp, "betad"), r$ri$entries)
  expect_length(sg, 1L)
  expect_identical(sg[[1L]]$dist$kind, "flat")
  nlp <- vapply(c(-0.8, 0, 0.9), function(a) {
    fp_nlp(fit, sg, list(betad = 0.3, thetaac = a))
  }, 0)
  expect_identical(nlp[1], nlp[2])
  expect_identical(nlp[1], nlp[3])
  lw <- r$ri$lower
  expect_length(lw, 1L)
  expect_lt(abs(unname(lw) - log(0.1)), 1e-12 * abs(log(0.1)))
})

test_that("a flat ar(cov = TRUE) run passes its prior_summary() back", {
  skip_sampler()
  d <- fp_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, prior = "flat", chains = 1, iter = 300, seed = 3,
               refresh = 0)))
  ps <- prior_summary(ds)
  e1 <- fp_rp(fit, "flat")$ri$entries
  set.seed(99)
  pts <- replicate(20, list(betad = stats::rnorm(1, 0, 0.7),
                            thetaac = stats::rnorm(1, 0, 0.7),
                            beta = stats::rnorm(2)), simplify = FALSE)
  for (back in list(ps, as.brmsprior(as.data.frame(ps)))) {
    e2 <- fp_rp(fit, back)$ri$entries
    v1 <- vapply(pts, function(p) fp_nlp(fit, e1, p), 0)
    v2 <- vapply(pts, function(p) fp_nlp(fit, e2, p), 0)
    expect_lt(max(abs(v1 - v2)), 1e-12 * max(abs(v1)))
  }
  # and frm()'s table of the same model is refused, as is its coerced
  # copy, which names the defaults its (flat) rows would replace
  p <- default_prior(y ~ x + ar(t, g, cov = TRUE), data = d)
  expect_error(fp_rp(fit, p), "table of route = \"fit\"", fixed = TRUE)
  expect_error(fp_rp(fit, as.brmsprior(as.data.frame(p))),
               "brms's default on Intercept, sigma with (flat)", fixed = TRUE)
})

test_that("(flat) on a mo() simplex is brms's dirichlet(1)", {
  # brms keeps dirichlet_lpdf(simo | con_simo) under "" with con_simo 1
  # (standata(): 1 1 1 on four values), so (flat) on the simplex is the
  # default. It was refused, "Invalid Dirichlet prior ... got 0"
  set.seed(3)
  dm <- data.frame(m2 = rep(1:2, 30), k = sample(1:4, 60, TRUE))
  dm$y <- 0.5 * dm$m2 + dm$k + stats::rnorm(60)
  fit <- frm(bf(y ~ mo(m2) + mo(k)), data = dm, dry_run = "objective")
  key <- function(es) sort(vapply(es, function(e) {
    paste(e$comp, paste(e$idx, collapse = ","), e$dist$kind)
  }, ""))
  e0 <- key(fp_rp(fit, NULL)$ri$entries)
  for (co in c("mom21", "mok1")) {
    e1 <- key(fp_rp(fit, set_prior("", class = "simo", coef = co))$ri$entries)
    expect_identical(e1, e0, label = co)
  }
  # the guard absent: a density of the right length still applies, and
  # one of the wrong length is still refused
  e2 <- key(fp_rp(fit, set_prior("dirichlet(2)", class = "simo",
                                  coef = "mok1"))$ri$entries)
  expect_true(any(grepl("dirichlet", e2)))
  expect_error(fp_rp(fit, set_prior("dirichlet(c(1, 2))", class = "simo",
                                    coef = "mok1")),
               "Invalid Dirichlet prior", fixed = TRUE)
})
