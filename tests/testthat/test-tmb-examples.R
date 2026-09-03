# The classic TMB examples, as frm() formulas.
#
# One row per example that frmtmb replicates. Each test builds the
# example's own negative log-likelihood inline with RTMB::MakeADFun(),
# optimizes it, and requires frm() to reach the same maximum. Both sides
# tape through RTMB, so the two objectives are the same arithmetic
# reached from two different directions: a disagreement is a modeling
# surface difference, not numerical noise.
#
# Nothing is downloaded and nothing is read from disk. Where an example
# ships a data file, the data is generated inline with the same
# structure and a fixed seed instead. Where an example is large
# (orange_big at 5000 latent variables, longlinreg at 10^6 rows) the
# generator is scaled down: the identity under test is the likelihood,
# not the row count. dev/tmb-examples-audit.md carries the full-size
# numbers and the verdict for every example in the upstream directory,
# including the ones no test here covers.
#
# @srrstats {G5.4a} Every row derives the answer twice: the example's
#   published likelihood, hand-written against RTMB, against frmtmb's
#   formula compiler.
skip_on_cran()

test_that("linreg: gaussian regression matches the TMB example", {
  set.seed(123)
  dd <- data.frame(Y = rnorm(10) + 1:10, x = 1:10)

  nll_ref <- function(p) {
    -sum(RTMB::dnorm(dd$Y, p$a + p$b * dd$x, exp(p$logSigma), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(a = 0, b = 0, logSigma = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  fit <- frm(bf(Y ~ x), family = gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  expect_vector_equal(unname(fixef(fit)$mu), unname(opt$par[c("a", "b")]),
                      tol = 1e-5)
})

test_that("dataeval: the pooled regression behind the tape-reuse demo", {
  # The example splits 100 rows into ten chunks and reuses one tape for
  # all of them, with a, b and sd shared. The fitted model is therefore
  # one pooled regression, which is what frm() writes.
  set.seed(124)
  dd <- data.frame(x = 1:100, y = 1:100 + rnorm(100))

  nll_ref <- function(p) {
    -sum(RTMB::dnorm(dd$y, p$a * dd$x + p$b, p$sd, log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(a = 0, b = 0, sd = 1), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = c(-Inf, -Inf, 1e-8))

  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("tweedie: the three-parameter Tweedie matches", {
  # The example reads a stored draw from tweedie::rtweedie. A compound
  # Poisson-gamma sum is the same distribution for 1 < p < 2 and needs
  # no extra package.
  set.seed(1001)
  n <- 1000
  mu <- 2; phi <- 2; pw <- 1.5
  lambda <- mu^(2 - pw) / (phi * (2 - pw))
  shape <- (2 - pw) / (pw - 1)
  scale <- phi * (pw - 1) * mu^(pw - 1)
  N <- rpois(n, lambda)
  y <- vapply(N, function(k) if (k == 0) 0 else sum(rgamma(k, shape, scale = scale)),
              numeric(1))
  dd <- data.frame(y = y)

  nll_ref <- function(p) {
    -sum(RTMB::dtweedie(dd$y, p$mu, p$phi, p$p, log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(mu = 1.1, phi = 1.1, p = 1.1),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = c(1e-8, 1e-8, 1 + 1e-8), upper = c(Inf, Inf, 2 - 1e-8))

  fit <- frm(bf(y ~ 1), family = tweedie(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("compois: the mean-parameterized Conway-Maxwell-Poisson matches", {
  # frmtmb's compois() is RTMB::dcompois2, the MEAN parameterization,
  # which is the example's second fit. The mode parameterization
  # (dcompois) has no family spelling.
  set.seed(123)
  nu <- .1; mode <- 10; domain <- 0:100
  prob <- dpois(domain, lambda = mode)^nu
  prob <- prob / sum(prob)
  x <- sample(domain, size = 2000, replace = TRUE, prob = prob)
  dd <- data.frame(x = x)

  nll_ref <- function(p) {
    -sum(RTMB::dcompois2(dd$x, exp(p$logmu), exp(p$lognu), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(logmu = 0, lognu = 0), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  fit <- frm(bf(x ~ 1), family = compois(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("spatial: a Poisson GLMM over an exponential covariance field", {
  # exp(pos + 0 | grp) is Sigma = sd^2 exp(-D / rho), which is the
  # example's exp(log_sigma) * u with u ~ MVN(0, exp(-a D)) and
  # rho = 1 / a. The example's own data: counts on a 10 x 10 lattice
  # with one covariate. The surface is flat enough that a resimulated
  # data set puts the two optimizers 3e-05 apart, so the covariate is
  # inlined verbatim and the positions are the lattice they are.
  sp_y <- c(2, 5, 0, 2, 4, 2, 8, 6, 8, 5, 0, 2, 2, 4, 0, 3, 9, 8, 4, 4, 1,
    3, 3, 3, 1, 2, 4, 5, 6, 2, 3, 3, 4, 4, 1, 5, 6, 4, 1, 1, 2, 1, 3, 11,
    2, 4, 13, 4, 4, 10, 3, 6, 1, 1, 0, 8, 4, 4, 5, 3, 1, 1, 1, 1, 4, 3, 1,
    4, 4, 4, 0, 0, 3, 7, 5, 4, 2, 2, 5, 2, 3, 1, 6, 6, 4, 3, 0, 5, 4, 3, 3,
    6, 4, 1, 5, 0, 5, 2, 5, 2)
  sp_x2 <- c(-0.07843052781, 0.04310351427, -0.1048777162, -0.06536390528,
    -0.002073643334, 0.1417429733, 0.1491268465, 0.02481202881,
    0.02395858034, 0.04416646603, 0.08421452044, 0.07816462156,
    -0.01402867919, 0.1080776553, 0.0693553828, -0.02497216501,
    -0.08544162512, -0.104265397, -0.1132628084, -0.1482812231,
    -0.1278271045, -0.0413755412, -0.09773557081, 0.1166788372,
    0.1450642589, 0.08466871813, -0.03733735931, 0.005612015261,
    -0.05610504542, 0.05180862542, 0.1209554956, 0.09795381097,
    -0.1169624834, -0.04232377996, 0.08903943555, 0.1359848434,
    -0.1441505658, -0.09405978005, 0.1242924539, 0.105890163,
    0.09817684095, -0.1418812039, 0.1226480866, 0.04612162634,
    0.0248917492, -0.02593747267, -0.01198069575, -0.1559568073,
    -0.1361770874, -0.06191070786, -0.07728992019, 0.03067909238,
    0.05207761075, -0.1420807121, -0.08397627334, -0.09647075431,
    0.1254553209, -0.1335876846, 0.0885559762, 0.1337777287, -0.136156121,
    -0.1365635438, 0.0888431759, 0.04177757569, -0.1500563082,
    -0.1470516612, 0.04021846004, 0.05352493282, 0.0356723245,
    -0.03024773017, 0.141994733, 0.05933785783, 0.00262609269,
    0.1052381875, -0.1592677825, 0.1024408573, -0.1437906496,
    0.1251007492, 0.1024546471, 0.05766970809, -0.1345753491,
    0.06676644743, -0.1535758012, 0.1380653958, 0.1376504082,
    -0.001600938292, -0.154250021, 0.1487791286, -0.1585675798,
    0.02025547526, 0.07217634359, 0.08140378473, -0.02090239182,
    -0.1115905231, -0.1259279076, 0.1188467067, 0.09626525633,
    0.005038780197, 0.01927707001, -0.09823082855)
  n <- 100
  Z <- cbind(rep(1:10, each = 10), rep(1:10, times = 10))
  dmat <- as.matrix(dist(Z))
  dd <- data.frame(y = sp_y, x2 = sp_x2,
                   pos = num_factor(Z[, 1], Z[, 2]),
                   grp = factor(rep(1L, n)))

  nll_ref <- function(p) {
    eta <- p$b[1] + p$b[2] * dd$x2 + exp(p$log_sigma) * p$u
    -RTMB::dmvnorm(p$u, 0, exp(-p$a * dmat), log = TRUE) -
      sum(RTMB::dpois(dd$y, exp(eta), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b = c(0, 0), a = 1.428571,
                              log_sigma = -0.6931472, u = rep(0, n)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = c(-100, -100, 0.01, -3),
                upper = c(100, 100, 3, 3))

  fit <- frm(bf(y ~ 1 + x2 + exp(pos + 0 | grp)), family = poisson(),
             data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("spde: Weibull survival with right censoring over an SPDE field", {
  skip_if_not_installed("fmesher")
  set.seed(126)
  n <- 200
  loc <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  mesh <- fmesher::fm_mesh_2d(loc = loc, max.edge = c(0.2, 0.5),
                              cutoff = 0.05)
  fem <- fmesher::fm_fem(mesh)
  fem <- list(c0 = fem$c0, g1 = fem$g1, g2 = fem$g2)
  idx <- mesh$idx$loc

  # Draw the field from the model's own precision, otherwise the SPDE
  # parameters are unidentified and the two optimizers wander apart on a
  # flat ridge instead of disagreeing about the likelihood.
  tau_t <- exp(-1); kap_t <- exp(1.5)
  Q_t <- tau_t^2 * (kap_t^4 * fem$c0 + 2 * kap_t^2 * fem$g1 + fem$g2)
  x_t <- as.vector(backsolve(chol(as.matrix(Q_t)), rnorm(nrow(fem$c0))))

  sexv <- rbinom(n, 1, 0.5)
  agev <- rnorm(n)
  eta_true <- -0.5 + 0.3 * sexv + 0.2 * agev + x_t[idx]
  time <- rweibull(n, shape = 1.2, scale = exp(-eta_true))
  cens <- as.integer(time > 2)          # brms convention: 1 = right censored
  time <- pmin(time, 2)
  dd <- data.frame(time = time, cens = cens, sex = sexv, age = agev,
                   node = factor(idx, levels = seq_len(mesh$n)))

  nll_ref <- function(p) {
    tau <- exp(p$log_tau); kap <- exp(p$log_kappa)
    Q <- tau^2 * (kap^4 * fem$c0 + 2 * kap^2 * fem$g1 + fem$g2)
    om <- exp(p$log_omega)
    eta <- p$beta[1] + p$beta[2] * dd$sex + p$beta[3] * dd$age +
      p$x[idx]
    lam <- exp(eta)
    # ifelse() would evaluate both branches and strip the advector
    # class; index the two groups instead, as the example's loop does
    ob <- which(dd$cens == 0)
    cn <- which(dd$cens == 1)
    ll_ob <- RTMB::dweibull(dd$time[ob], shape = om, scale = 1 / lam[ob],
                            log = TRUE)
    ll_cn <- log(1 - RTMB::pweibull(dd$time[cn], shape = om,
                                    scale = 1 / lam[cn]))
    -RTMB::dgmrf(p$x, 0, Q, log = TRUE) - sum(ll_ob) - sum(ll_cn)
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(beta = c(-0.5, 0, 0), log_tau = -2,
                              log_kappa = 1, log_omega = 0,
                              x = rep(0, nrow(fem$c0))),
                         random = "x", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  fit <- frm(bf(time | cens(cens) ~ sex + age + spde(fem, gr = node)),
             family = weibull(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-4)
})

test_that("adaptive_integration: exact per-observation marginalization", {
  # Both sides marginalize each scalar random effect exactly. The
  # reference uses stats::integrate (adaptive QUADPACK over the line);
  # frm(quadrature = TRUE) uses TMB's fixed Gauss-Kronrod rule. The
  # tolerance below is the rule difference, not a model difference: it
  # is 7.9e-07 of the objective. See dev/tmb-examples-audit.md.
  set.seed(123)
  ndat <- 300
  c1 <- rnorm(ndat); c2 <- rnorm(ndat)
  ngroup <- 5
  c3 <- cut(runif(ndat), ngroup)
  eps <- rnorm(ndat, sd = .5)
  pr <- plogis(c1 + c2 + seq(-1, 1, length = ngroup)[c3] + eps)
  nn <- rep(10, ndat)
  xx <- rbinom(ndat, size = nn, prob = pr)
  A <- model.matrix(~ c1 + c2 + c3 - 1)
  dd <- data.frame(x = xx, n = nn, c1 = c1, c2 = c2, c3 = c3,
                   obs = factor(seq_len(ndat)))

  # Vectorize() over advectors needs RTMB attached to simplify its
  # mapply result without losing the class; do.call("c", lapply(...))
  # dispatches on the registered method and needs no search-path change.
  gauss_binom <- function(x, n, mu, sd) {
    RTMB::integrate(function(u) {
      exp(RTMB::dnorm(u, 0, 1, log = TRUE) +
            RTMB::dbinom_robust(x, n, sd * u + mu, log = TRUE))
    }, -Inf, Inf)$value
  }
  nll_ref <- function(p) {
    mu <- as.vector(A %*% p$b)
    sd <- exp(p$logsd)
    v <- do.call("c", lapply(seq_len(ndat), function(i) {
      gauss_binom(dd$x[i], dd$n[i], mu[i], sd)
    }))
    -sum(log(v))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = rep(0, ncol(A)), logsd = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  fit <- frm(bf(x | trials(n) ~ 0 + c1 + c2 + c3 + (1 | obs)),
             family = binomial(), data = dd, quadrature = TRUE)
  rel <- abs(as.numeric(logLik(fit)) - (-opt$objective)) / abs(opt$objective)
  expect_lt(rel, 1e-5)
})

test_that("transform: a latent AR(1) field pushed through a gamma quantile", {
  # dautoreg() fixes the field's scale at 1, because pnorm(u) is only a
  # uniform when u is marginally standard normal. frmtmb's ar1() block
  # carries a free marginal sd, so theta_1 is pinned to log(1) = 0 with
  # lower/upper. That is the map= of the reference, spelled in bounds.
  set.seed(123)
  n <- 200; phi <- .6
  u <- numeric(n); u[1] <- rnorm(1)
  for (i in 2:n) u[i] <- phi * u[i - 1] + rnorm(1, sd = sqrt(1 - phi^2))
  y <- qgamma(pnorm(u), shape = 2, scale = 3) + rnorm(n, sd = 2)
  dd <- data.frame(y = y, tim = factor(seq_len(n)),
                   g = factor(rep(1L, n)))

  nll_ref <- function(p) {
    -RTMB::dautoreg(p$u, phi = p$phi, log = TRUE) -
      sum(RTMB::dnorm(y, RTMB::qgamma(RTMB::pnorm(p$u), p$shape,
                                      scale = p$scale),
                      p$sd, log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(phi = 0, shape = 1, scale = 1, sd = 1,
                              u = u * 0),
                         random = "u", silent = TRUE)
  opt <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr, lower = 1e-6))

  # The nl body is evaluated in the formula's own environment, so the
  # RTMB:: prefixes are load bearing: bare qgamma/pnorm find the stats
  # versions unless the user has attached RTMB.
  pin <- c(theta_1 = 0)
  fit <- suppressWarnings(
    frm(bf(y ~ RTMB::qgamma(RTMB::pnorm(z), shape, scale),
           z ~ 0 + ar1(tim + 0 | g), shape ~ 1, scale ~ 1, nl = TRUE),
        family = gaussian(), data = dd, start = list(beta = c(0, 2, 3)),
        lower = pin, upper = pin))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  # theta_2 -> rho is the reference's phi
  th <- fit$estimates$theta
  expect_lt(abs(th[2] / sqrt(1 + th[2]^2) - opt$par[["phi"]]), 1e-4)
})

test_that("transform2: the same field through a beta quantile", {
  set.seed(123)
  n <- 200; phi <- .6
  u <- numeric(n); u[1] <- rnorm(1)
  for (i in 2:n) u[i] <- phi * u[i - 1] + rnorm(1, sd = sqrt(1 - phi^2))
  y <- qbeta(pnorm(u), shape1 = .5, shape2 = 2) + rnorm(n, sd = .005)
  dd <- data.frame(y = y, tim = factor(seq_len(n)),
                   g = factor(rep(1L, n)))

  nll_ref <- function(p) {
    -RTMB::dautoreg(p$u, phi = p$phi, log = TRUE) -
      sum(RTMB::dnorm(y, RTMB::qbeta(RTMB::pnorm(p$u), p$shape1, p$shape2),
                      p$sd, log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(phi = 0, shape1 = 1, shape2 = 1, sd = 1,
                              u = u * 0),
                         random = "u", silent = TRUE)
  opt <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr, lower = 1e-6))

  pin <- c(theta_1 = 0)
  fit <- suppressWarnings(
    frm(bf(y ~ RTMB::qbeta(RTMB::pnorm(z), shape1, shape2),
           z ~ 0 + ar1(tim + 0 | g), shape1 ~ 1, shape2 ~ 1, nl = TRUE),
        family = gaussian(), data = dd, start = list(beta = c(0, .5, 2)),
        lower = pin, upper = pin))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("orange_big: logistic growth with a random asymptote", {
  # The example's data, replicated `mult` times to make one latent
  # random effect per replicated tree. The upstream run uses
  # mult = 1000 (5000 latent variables); the audit reports that size,
  # the test uses a tenth of it.
  yv <- c(30, 58, 87, 115, 120, 142, 145, 33, 69, 111, 156, 172, 203, 203,
          30, 51, 75, 108, 115, 139, 140, 32, 62, 112, 167, 179, 209, 214,
          30, 49, 81, 125, 142, 174, 177)
  tv <- rep(c(118, 484, 664, 1004, 1231, 1372, 1582), 5)
  mult <- 100
  dd <- data.frame(y = rep(yv, mult), t = rep(tv, mult),
                   tree = factor(rep(seq_len(5 * mult), each = 7)))

  nll_ref <- function(p) {
    a0 <- 192 + p$beta[1] + p$u[as.integer(dd$tree)]
    f <- a0 / (1 + exp(-(dd$t - (726 + p$beta[2])) / (356 + p$beta[3])))
    -sum(RTMB::dnorm(p$u, 0, exp(p$log_sigma_u), log = TRUE)) -
      sum(RTMB::dnorm(dd$y, f, exp(p$log_sigma), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(beta = c(0, 0, 0), log_sigma = 1,
                              log_sigma_u = 2, u = rep(0, 5 * mult)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = c(-10, -10, -10, -5, -5),
                upper = c(10, 10, 10, 5, 5))

  fit <- suppressWarnings(
    frm(bf(y ~ a0 / (1 + exp(-(t - a1) / a2)),
           a0 ~ 1 + (1 | tree), a1 ~ 1, a2 ~ 1, nl = TRUE),
        family = gaussian(), data = dd,
        start = list(beta = c(192, 726, 356))))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-5)
  expect_lt(abs(fixef(fit)$a0[[1]] - (192 + opt$par[1])), 1e-3)
})

test_that("socatt: cumulative logit with a random intercept", {
  # The ADMB original writes the group effect as sigma * u with
  # u ~ N(0, 1); (1 | g) is the same block in its own parameterization,
  # so no pin is needed.
  set.seed(127)
  M <- 200                       # groups, 4 observations each
  S <- 5                         # response categories
  g <- factor(rep(seq_len(M), each = 4))
  n <- length(g)
  X <- cbind(rbinom(n, 1, .4), rnorm(n), runif(n))
  alpha_true <- c(-1.5, -0.4, 0.5, 1.6)
  eta <- as.vector(X %*% c(0.6, -0.3, 0.4)) + rnorm(M, 0, 0.7)[g]
  cum <- cbind(vapply(alpha_true, function(a) plogis(a - eta), numeric(n)), 1)
  yv <- vapply(seq_len(n), function(i) findInterval(runif(1), cum[i, ]) + 1L,
               integer(1))
  yv <- pmin(pmax(yv, 1L), S)
  dd <- data.frame(y = ordered(yv), g = g,
                   x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    alpha <- p$tmpk
    for (s in 2:length(p$tmpk)) alpha[s] <- alpha[s - 1] + exp(p$tmpk[s])
    et <- as.vector(X %*% p$b) + exp(p$logsigma) * p$u[as.integer(g)]
    P <- rep(0, n)
    hi <- yv < S
    P[hi] <- 1 / (1 + exp(-(alpha[yv[hi]] - et[hi])))
    P[!hi] <- 1
    lo <- yv > 1
    P[lo] <- P[lo] - 1 / (1 + exp(-(alpha[yv[lo] - 1] - et[lo])))
    -sum(RTMB::dnorm(p$u, 0, 1, log = TRUE)) - sum(log(P))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b = rep(0, 3), logsigma = -0.5,
                              tmpk = rep(0, S - 1), u = rep(0, M)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  fit <- frm(bf(y ~ x1 + x2 + x3 + (1 | g)), family = cumulative(),
             data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-5)
})

test_that("lr_test: the map= restrictions are formula restrictions", {
  # The example uses map= to fit three nested models on a ragged array
  # and runs two likelihood-ratio tests. In frmtmb the restriction is a
  # formula, and the parameter counts agree, so anova() reproduces the
  # example's table.
  ngroup <- 5
  nrep <- c(5, 8, 11, 13, 2)
  sd_true <- c(1, 1, 1, 2, 2)
  set.seed(123)
  ra <- lapply(seq_len(ngroup), function(i) rnorm(nrep[i], 0, sd_true[i]))
  dd <- data.frame(obs = unlist(ra),
                   g = factor(rep(seq_len(ngroup), lengths(ra))))

  nll_ref <- function(p) {
    gi <- as.integer(dd$g)
    -sum(RTMB::dnorm(dd$obs, p$mu[gi], p$sd[gi], log = TRUE))
  }
  ref_ll <- function(map) {
    o <- RTMB::MakeADFun(nll_ref, list(mu = rep(0, 5), sd = rep(1, 5)),
                         map = map, silent = TRUE)
    op <- suppressWarnings(nlminb(o$par, o$fn, o$gr))
    list(ll = -op$objective, npar = length(o$par))
  }
  r_full <- ref_ll(NULL)
  r_1 <- ref_ll(list(mu = factor(rep(1, 5))))
  r_2 <- ref_ll(list(mu = factor(rep(1, 5)), sd = factor(rep(1, 5))))

  f_full <- frm(bf(obs ~ 0 + g, sigma ~ 0 + g), family = gaussian(),
                data = dd)
  f_1 <- frm(bf(obs ~ 1, sigma ~ 0 + g), family = gaussian(), data = dd)
  f_2 <- frm(bf(obs ~ 1, sigma ~ 1), family = gaussian(), data = dd)

  expect_lt(abs(as.numeric(logLik(f_full)) - r_full$ll), 1e-6)
  expect_lt(abs(as.numeric(logLik(f_1)) - r_1$ll), 1e-6)
  expect_lt(abs(as.numeric(logLik(f_2)) - r_2$ll), 1e-6)
  expect_equal(length(f_full$opt$par), r_full$npar)
  expect_equal(length(f_1$opt$par), r_1$npar)
  expect_equal(length(f_2$opt$par), r_2$npar)
})

test_that("longlinreg: the same regression at scale", {
  # The upstream example is 10^6 rows; the audit reports that size
  # (7.9 s reference against 8.4 s frm(), agreeing to 1.86e-09). The
  # test runs a tenth of it.
  set.seed(123)
  nobs <- 1e5
  x <- seq(0, 10, length = nobs)
  dd <- data.frame(Y = 2 * x + 1 + rnorm(nobs), x = x)

  nll_ref <- function(p) {
    -sum(RTMB::dnorm(dd$Y, p$a + p$b * dd$x, exp(p$logSigma), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(a = 0, b = 0, logSigma = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)

  fit <- suppressWarnings(frm(bf(Y ~ x), family = gaussian(), data = dd))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)) /
              abs(opt$objective), 1e-10)
})

test_that("sdv_multi: multivariate stochastic volatility, at a short length", {
  # AWKWARD in the audit, but only for cost: the grammar maps exactly.
  # rescor gives the example's correlation matrix R, sigma's log link
  # makes its intercept mu_x / 2 and its latent block h / 2, and h / 2
  # is AR(1) with the same phi, so ar1()'s marginal sd absorbs the
  # factor of 2. What does not scale is ar1() itself: it builds a dense
  # d x d covariance, so the example's own d = 945 does not tape. This
  # runs at d = 40, where the identity is exact.
  p <- 3
  # sdv_multi.R's explicit Cholesky parameterization of the correlation
  # matrix, rather than RTMB::unstructured(), whose helper needs RTMB on
  # the search path to keep the advector class through its sub-assignment
  corr_of <- function(od) {
    L <- RTMB::matrix(c(1, od[1], od[2], 0, 1, od[3], 0, 0, 1), 3, 3)
    rn <- sqrt(c(1, od[1]^2 + 1, od[2]^2 + od[3]^2 + 1))
    L <- L / rn
    L %*% t(L)
  }
  set.seed(77)
  n <- 40
  phi_t <- c(0.9, 0.85, 0.8); sig_t <- c(0.3, 0.35, 0.25)
  mux_t <- c(-0.5, -0.4, -0.6)
  h <- matrix(0, n, p)
  for (j in seq_len(p)) {
    h[1, j] <- rnorm(1, 0, sig_t[j] / sqrt(1 - phi_t[j]^2))
    for (i in 2:n) h[i, j] <- phi_t[j] * h[i - 1, j] + rnorm(1, 0, sig_t[j])
  }
  Rc <- matrix(c(1, .4, .2, .4, 1, .3, .2, .3, 1), 3)
  Xd <- matrix(0, n, p)
  for (i in seq_len(n)) {
    sy <- exp(0.5 * (mux_t + h[i, ]))
    Xd[i, ] <- as.vector(t(chol(diag(sy) %*% Rc %*% diag(sy))) %*% rnorm(p))
  }

  nll_ref <- function(pp) {
    hh <- t(pp$h)
    si <- exp(pp$log_sigma) / sqrt(1 - pp$phi^2)
    nll <- 0
    for (j in seq_len(p)) {
      nll <- nll - RTMB::dautoreg(hh[, j], phi = pp$phi[j], scale = si[j],
                                  log = TRUE)
    }
    # base matrix() strips the advector class; RTMB::matrix keeps it
    mm <- RTMB::matrix(pp$mu_x, nrow(hh), ncol(hh), byrow = TRUE)
    nll - sum(RTMB::dmvnorm(Xd, 0, corr_of(pp$off_diag_x),
                            scale = exp(.5 * (hh + mm)), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(phi = rep(0.9, p), log_sigma = rep(-1.7, p),
                              mu_x = rep(-0.5, p),
                              off_diag_x = rep(0, 3),
                              h = matrix(0, p, n)),
                         random = "h", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = c(rep(-.99, 3), rep(-3, 9)),
                upper = c(rep(.99, 3), rep(3, 9)))

  dd <- data.frame(x1 = Xd[, 1], x2 = Xd[, 2], x3 = Xd[, 3],
                   tim = factor(seq_len(n)), g = factor(rep(1L, n)))
  fit <- suppressWarnings(
    frm(mvbf(bf(x1 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
             bf(x2 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
             bf(x3 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
             rescor = TRUE), data = dd))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("mvrw: Laplace and a Kalman filter marginalize the same model", {
  # NOTE: this block asserts an RTMB-vs-RTMB identity (Laplace vs the
  # Kalman filter on the same model) and calls no frmtmb code. It is
  # deliberate dependency insurance for the roadmap decision that no
  # Kalman node is needed; a failure here means upstream RTMB moved,
  # not frmtmb.

  # mvrw is out of scope for frm() (no latent random-walk term), but the
  # identity behind the Kalman-node prototype is cheap to hold onto: a
  # linear-Gaussian state space marginalized by Laplace over the states
  # and by a Kalman filter in closed form must agree at every parameter
  # value. dev/tmbex-kalman-prototype.R has the tape sizes and timings.
  skip_if_not_installed("MASS")
  set.seed(1)
  sd_dim <- 3; nT <- 60
  Sig <- 0.9^abs(outer(1:3, 1:3, "-")) *
    (seq(0.5, 2, length = 3) %o% seq(0.5, 2, length = 3))
  L <- t(chol(Sig))
  st <- matrix(0, sd_dim, nT)
  st[, 1] <- rnorm(sd_dim)
  for (i in 2:nT) st[, i] <- st[, i - 1] + L %*% rnorm(sd_dim)
  Y <- st + matrix(rnorm(sd_dim * nT), sd_dim)

  trf <- function(x) 2 / (1 + exp(-2 * x)) - 1
  pieces <- function(th) {
    rho <- trf(th[1]); sds <- exp(th[2:4]); so <- exp(th[5:7])
    list(Sigma = outer(1:3, 1:3,
                       function(i, j) rho^abs(i - j) * sds[i] * sds[j]),
         R = so^2 * diag(3))
  }
  joint <- function(p) {
    pc <- pieces(c(p$transf_rho, p$logsds, p$logsdObs))
    -sum(RTMB::dmvnorm(diff(t(p$u)), 0, pc$Sigma, log = TRUE)) -
      sum(RTMB::dnorm(Y, p$u, exp(p$logsdObs), log = TRUE))
  }
  # Flat prior on the first state: integrating N(y1; u1, R) over u1 is
  # exactly 1, so the filter starts at a = y1, P = R and observation one
  # contributes nothing. That is what makes the two routes agree
  # constant and all.
  kalman <- function(th) {
    pc <- pieces(th)
    a <- Y[, 1]; P <- pc$R; nll <- 0
    for (t in 2:nT) {
      P <- P + pc$Sigma
      v <- Y[, t] - a
      Fm <- P + pc$R
      nll <- nll - RTMB::dmvnorm(v, 0, Fm, log = TRUE)
      K <- P %*% RTMB::solve(Fm)
      a <- a + as.vector(K %*% v)
      P <- P - K %*% P
    }
    nll
  }
  th0 <- c(0.1, rep(0, 6))
  obj_l <- RTMB::MakeADFun(joint,
                           list(transf_rho = 0.1, logsds = rep(0, 3),
                                logsdObs = rep(0, 3), u = Y * 0),
                           random = "u", silent = TRUE)
  TMB::newtonOption(obj_l, smartsearch = FALSE)
  obj_k <- RTMB::MakeADFun(function(p) kalman(p$th), list(th = th0),
                           silent = TRUE)

  # the identity, away from the optimum as well as at it
  for (th in list(th0,
                  c(0.5, 0.2, -0.3, 0.4, -0.1, 0.2, 0.0),
                  c(-0.4, -0.5, 0.6, 0.1, 0.3, -0.2, 0.5))) {
    expect_lt(abs(obj_l$fn(th) - obj_k$fn(th)), 1e-8)
  }
  opt_l <- nlminb(obj_l$par, obj_l$fn, obj_l$gr)
  opt_k <- nlminb(obj_k$par, obj_k$fn, obj_k$gr)
  expect_lt(abs(opt_l$objective - opt_k$objective), 1e-8)
  expect_lt(max(abs(opt_l$par - opt_k$par)), 1e-5)
})
