# data2 = list(): structural objects (adjacency matrices, precisions,
# covariances, mesh triples) handed to the fit by value instead of being
# captured from the calling environment.

# Rook-adjacency of an r x c lattice, with brms-style dimnames.
d2_lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c))
  n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
    }
  }
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}

# alpha = 2 finite-element triple of a regular 1-D linear mesh.
d2_chain_fem <- function(nn, h = 0.5) {
  C0 <- diag(rep(h, nn))
  C0[1, 1] <- h / 2
  C0[nn, nn] <- h / 2
  G <- matrix(0, nn, nn)
  for (i in seq_len(nn - 1)) {
    G[i, i] <- G[i, i] + 1 / h
    G[i + 1, i + 1] <- G[i + 1, i + 1] + 1 / h
    G[i, i + 1] <- -1 / h
    G[i + 1, i] <- -1 / h
  }
  list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
}

# A tri-diagonal precision over ng levels plus a data set drawn from it.
d2_prec_data <- function(seed, ng = 8, per = 10) {
  set.seed(seed)
  Q <- Matrix::bandSparse(ng, k = c(-1, 0, 1),
                          diagonals = list(rep(-0.4, ng - 1),
                                           rep(1.2, ng),
                                           rep(-0.4, ng - 1)))
  dimnames(Q) <- list(paste0("g", seq_len(ng)), paste0("g", seq_len(ng)))
  A <- solve(as.matrix(Q))
  dimnames(A) <- dimnames(Q)
  b <- drop(crossprod(chol(A), stats::rnorm(ng))) * 0.7
  dd <- data.frame(g = factor(rep(rownames(A), each = per),
                              levels = rownames(A)),
                   x = stats::rnorm(ng * per))
  dd$y <- 1 + 0.5 * dd$x + b[as.integer(dd$g)] +
    stats::rnorm(nrow(dd), 0, 0.5)
  list(dd = dd, Q = Q, A = A)
}

# The data2 spelling always uses a name that exists NOWHERE else, so a
# fit that resolved it could only have read data2.
test_that("every structural special resolves from data2", {
  s <- d2_prec_data(11)
  dd <- s$dd
  Q <- s$Q
  A <- s$A

  f_prec <- frm(bf(y ~ x + (1 | gr(g, prec = Q))) + gaussian(), data = dd)
  f_prec2 <- frm(bf(y ~ x + (1 | gr(g, prec = Qd2))) + gaussian(),
                 data = dd, data2 = list(Qd2 = Q))
  expect_equal(as.numeric(logLik(f_prec2)), as.numeric(logLik(f_prec)),
               tolerance = 1e-10)

  f_cov <- frm(bf(y ~ x + (1 | gr(g, cov = A))) + gaussian(), data = dd)
  f_cov2 <- frm(bf(y ~ x + (1 | gr(g, cov = Ad2))) + gaussian(),
                data = dd, data2 = list(Ad2 = A))
  expect_equal(as.numeric(logLik(f_cov2)), as.numeric(logLik(f_cov)),
               tolerance = 1e-10)

  set.seed(93)
  ng <- 30
  V <- matrix(c(1.2, 0.5, 0.5, 0.8), 2, 2)
  bb <- t(chol(V)) %*% matrix(stats::rnorm(2 * ng), 2)
  ed <- data.frame(y = 1 + as.vector(bb) + stats::rnorm(2 * ng, 0, 0.5),
                   f = factor(rep(c("a", "b"), ng)),
                   g = factor(rep(seq_len(ng), each = 2)))
  f_eq <- frm(bf(y ~ 1 + equalto(f + 0 | g, V)) + gaussian(), data = ed)
  f_eq2 <- frm(bf(y ~ 1 + equalto(f + 0 | g, Vd2)) + gaussian(),
               data = ed, data2 = list(Vd2 = V))
  expect_equal(as.numeric(logLik(f_eq2)), as.numeric(logLik(f_eq)),
               tolerance = 1e-10)

  set.seed(12)
  W <- d2_lattice_W(2, 3)
  nl <- nrow(W)
  cd <- data.frame(loc = factor(rep(rownames(W), each = 6),
                                levels = rownames(W)),
                   x = stats::rnorm(nl * 6))
  cd$y <- 1 + 0.4 * cd$x + stats::rnorm(nl, 0, 0.8)[as.integer(cd$loc)] +
    stats::rnorm(nrow(cd), 0, 0.5)
  f_car <- frm(bf(y ~ x + car(W, gr = loc)) + gaussian(), data = cd)
  f_car2 <- frm(bf(y ~ x + car(Wd2, gr = loc)) + gaussian(), data = cd,
                data2 = list(Wd2 = W))
  expect_equal(as.numeric(logLik(f_car2)), as.numeric(logLik(f_car)),
               tolerance = 1e-10)

  set.seed(13)
  nn <- 8
  fem <- d2_chain_fem(nn)
  sd_ <- data.frame(node = factor(rep(seq_len(nn), each = 6)),
                    x = stats::rnorm(nn * 6))
  sd_$y <- 0.3 + 0.4 * sd_$x +
    stats::rnorm(nn, 0, 0.7)[as.integer(sd_$node)] +
    stats::rnorm(nrow(sd_), 0, 0.3)
  f_spde <- frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(), data = sd_)
  f_spde2 <- frm(bf(y ~ x + spde(femd2, gr = node)) + gaussian(),
                 data = sd_, data2 = list(femd2 = fem))
  expect_equal(as.numeric(logLik(f_spde2)), as.numeric(logLik(f_spde)),
               tolerance = 1e-10)
})

test_that("data2 resolves compound expressions, unlike brms", {
  s <- d2_prec_data(14, ng = 6, per = 8)
  dd <- s$dd
  Q <- as.matrix(s$Q)
  f_cov <- frm(bf(y ~ x + (1 | gr(g, cov = A))) + gaussian(),
               data = dd, data2 = list(A = s$A))
  # brms takes a bare name from data2 and nothing else; here the whole
  # expression is evaluated with data2 in front of the data mask
  f_expr <- frm(bf(y ~ x + (1 | gr(g, cov = solve(Qd2)))) + gaussian(),
                data = dd, data2 = list(Qd2 = Q))
  expect_equal(as.numeric(logLik(f_expr)), as.numeric(logLik(f_cov)),
               tolerance = 1e-8)
})

test_that("data2 shadows a same-named column of data", {
  s <- d2_prec_data(15, ng = 6, per = 8)
  dd <- s$dd
  # a data column named A would win the old data-then-env lookup; data2
  # is searched first, so the matrix is what reaches gr(cov = )
  dd$A <- stats::rnorm(nrow(dd))
  f <- frm(bf(y ~ x + (1 | gr(g, cov = A))) + gaussian(), data = dd,
           data2 = list(A = s$A))
  expect_equal(f$frame$re_blocks[[1]]$n_levels, 6L)
})

test_that("a structural object found nowhere names data2", {
  s <- d2_prec_data(16, ng = 5, per = 6)
  expect_error(
    frm(bf(y ~ x + (1 | gr(g, prec = Qnope))) + gaussian(), data = s$dd,
        data2 = list(Q = s$Q)),
    "cannot find 'Qnope'"
  )
  expect_error(
    frm(bf(y ~ x + (1 | gr(g, prec = Qnope))) + gaussian(), data = s$dd,
        data2 = list(Q = s$Q)),
    "data2 = list\\(Qnope = Qnope\\)"
  )
  expect_error(
    frm(bf(y ~ x + (1 | gr(g, prec = Q))) + gaussian(), data = s$dd,
        data2 = list(s$Q)),
    "named list"
  )
  expect_error(
    frm(bf(y ~ x + (1 | gr(g, prec = Q))) + gaussian(), data = s$dd,
        data2 = c(1, 2)),
    "named list"
  )
})

test_that("data2 outlives its calling environment across saveRDS", {
  s <- d2_prec_data(17, ng = 6, per = 8)
  dd <- s$dd
  # the matrix lives only here, and the formula carries this
  # environment: it is the environment a pre-data2 fit would depend on
  env <- new.env(parent = globalenv())
  env$A <- s$A
  form <- y ~ x + (1 | gr(g, cov = A))
  environment(form) <- env

  fit_env <- frm(bf(form) + gaussian(), data = dd)
  fit_d2 <- frm(bf(form) + gaussian(), data = dd, data2 = list(A = s$A))
  expect_equal(as.numeric(logLik(fit_d2)), as.numeric(logLik(fit_env)),
               tolerance = 1e-10)

  # a fresh session has neither the binding nor the workspace that held
  # it; deleting it before serializing is what a new session looks like
  # from the lookup's side (the captured environment would otherwise be
  # written into the .rds along with the matrix it holds)
  rm("A", envir = env)
  f_env <- tempfile(fileext = ".rds")
  f_d2 <- tempfile(fileext = ".rds")
  saveRDS(fit_env, f_env)
  saveRDS(fit_d2, f_d2)
  rm(fit_env, fit_d2, env, form)
  gc(verbose = FALSE)

  r_env <- readRDS(f_env)
  r_d2 <- readRDS(f_d2)
  expect_equal(names(r_d2$data2), "A")

  # the data2 fit re-assembles; the environment-capture fit cannot
  expect_error(
    frmtmb:::assemble_frame(r_env$spec, dd),
    "cannot find 'A'"
  )
  expect_silent(frmtmb:::assemble_frame(r_d2$spec, dd,
                                        data2 = r_d2$data2))

  # influence() re-assembles once per deleted group: every refit of the
  # environment fit fails (NA rows), every refit of the data2 fit works
  i_d2 <- suppressWarnings(influence(r_d2, groups = "g"))
  expect_false(anyNA(i_d2$fixed))
  i_env <- suppressWarnings(influence(r_env, groups = "g"))
  expect_true(all(is.na(i_env$fixed)))

  # refit() reuses the assembled frame, so it only has to survive the
  # round trip: the restored objective is a dead pointer and is rebuilt
  rf <- refit(r_d2, stats::simulate(r_d2, nsim = 1, re.form = NA)[[1]])
  expect_s3_class(rf, "frmtmb_fit")
  expect_equal(names(rf$data2), "A")
  unlink(c(f_env, f_d2))
})

test_that("update() and drop1() re-assemble from the stored data2", {
  s <- d2_prec_data(18, ng = 6, per = 8)
  dd <- s$dd
  env <- new.env(parent = globalenv())
  env$A <- s$A
  form <- y ~ x + (1 | gr(g, cov = A))
  environment(form) <- env
  fit <- frm(bf(form) + gaussian(), data = dd, data2 = list(A = s$A))
  rm("A", envir = env)

  up <- stats::update(fit, data = dd[1:40, ])
  expect_equal(nobs(up), 40L)
  expect_equal(names(up$data2), "A")

  d1 <- stats::drop1(fit)
  expect_equal(rownames(d1), c("<none>", "x"))
  expect_false(anyNA(d1$AIC))
})

test_that("frm_multiple carries data2 into every imputation", {
  s <- d2_prec_data(19, ng = 6, per = 8)
  imps <- lapply(1:3, function(i) {
    d <- s$dd
    d$x <- d$x + stats::rnorm(nrow(d), 0, 0.05)
    d
  })
  env <- new.env(parent = globalenv())
  form <- y ~ x + (1 | gr(g, cov = A))
  environment(form) <- env
  m <- frm_multiple(bf(form) + gaussian(), data = imps,
                    data2 = list(A = s$A))
  expect_length(m$fits, 3L)
  for (f in m$fits) expect_equal(names(f$data2), "A")
  expect_true(all(is.finite(m$pooled$estimate)))
})

test_that("get_prior() and frm_simulate() take data2", {
  s <- d2_prec_data(20, ng = 6, per = 8)
  env <- new.env(parent = globalenv())
  form <- y ~ x + (1 | gr(g, cov = A))
  environment(form) <- env
  gp <- get_prior(bf(form) + gaussian(), data = s$dd,
                  data2 = list(A = s$A))
  expect_true("sd" %in% gp$class)
  sims <- frm_simulate(bf(form) + gaussian(), s$dd,
                       newparams = list(Intercept = 1, x = 0.5,
                                        sigma = 0.7,
                                        sd_g__Intercept = 0.5),
                       nsim = 2, seed = 1, data2 = list(A = s$A))
  expect_equal(dim(sims), c(nrow(s$dd), 2L))
})
