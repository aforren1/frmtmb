# Fits that stopped short of the maximum likelihood without saying so,
# and the optimizer paths that died on the way there (lane optima,
# dev/optima-findings.md).

# ---- mo(): the simplex chart and the search over the other sign -------

# brms_monotonic's own data code
optima_mo_data <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}

# The exact maximum of `ls ~ mo(income) * age`. For each sign of the two
# scale coefficients the increments b * D * w share that sign, so the
# model is least squares in them under sign constraints. The constrained
# optimum is the plain least-squares fit on its own support, so the best
# fit over every sign and support whose coefficients keep their signs is
# the maximum, exactly, with no start and no tolerance.
optima_mo_exact <- function(d) {
  code <- as.integer(d$income) - 1L
  S <- sapply(1:3, function(k) as.numeric(code >= k))
  n <- nrow(d)
  best <- -Inf
  for (s1 in c(-1, 1)) for (s2 in c(-1, 1)) for (m1 in 0:7) for (m2 in 0:7) {
    a1 <- which(bitwAnd(m1, c(1L, 2L, 4L)) > 0)
    a2 <- which(bitwAnd(m2, c(1L, 2L, 4L)) > 0)
    X <- cbind(1, d$age, S[, a1, drop = FALSE],
               (S * d$age)[, a2, drop = FALSE])
    f <- stats::lm.fit(X, d$ls)
    cf <- f$coefficients[-(1:2)]
    if (anyNA(cf)) next
    sg <- c(rep(s1, length(a1)), rep(s2, length(a2)))
    if (any(sg * cf < 0)) next
    best <- max(best,
                -(n / 2) * (log(2 * pi * sum(f$residuals^2) / n) + 1))
  }
  best
}

test_that("a mo() fit reaches the exact maximum", {
  # on 0.68.1, with code 0 and no warning: seed 12 stopped 1.95 below on
  # a softmax plateau, seed 54 1.58 below at the other sign's maximum,
  # and seed 38 0.28 below; on the simplex chart alone seed 38 crept
  # toward the chart's pole and stopped 0.23 below
  for (s in c(12, 54, 38)) {
    d <- optima_mo_data(s)
    fit <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
    ex <- optima_mo_exact(d)
    expect_identical(fit$opt$convergence, 0L)
    expect_lt((ex - as.numeric(logLik(fit))) / abs(ex), 1e-8,
              label = paste("relative gap below the maximum, seed", s))
  }
})

test_that("the simplex chart reaches every face at a finite point", {
  ms <- frmtmb:::mo_simplex
  mc <- frmtmb:::mo_coords
  expect_equal(ms(c(0, 0, 0)), rep(0.25, 4))
  set.seed(3)
  for (D in 2:5) {
    for (k in seq_len(D)) {
      v <- replace(numeric(D), k, 1)
      x <- mc(v)
      # a vertex, at a finite coordinate inside the unit ball
      expect_true(all(is.finite(x)))
      expect_lt(sum(x^2), 1)
      expect_equal(ms(x), v)
    }
    w <- stats::rexp(D)
    w <- w / sum(w)
    w[1] <- 0
    w <- w / sum(w)
    expect_equal(ms(mc(w)), w)
    expect_lt(sum(mc(w)^2), 1)
  }
  # the same simplex from the far side of the chart maps back inside
  z <- c(-4601.2, 3066.4)
  expect_gt(sum(z^2), 1)
  expect_equal(ms(mc(ms(z))), ms(z))
  expect_lt(sum(mc(ms(z))^2), 1)
})

test_that("an interior simplex is not searched, a face is", {
  # mo_terms_data(): the main effect's steps 0.7, 0.2, 0.1 are interior
  fit <- frm(bf(y ~ mo(inc) + z) + gaussian(), data = mo_terms_data())
  expect_null(fit$opt$mo_search)
  # seed 54's interaction sits on a face, so the other sign is tried.
  # What the search gains depends on where the chart alone stopped,
  # which moves with the rounding of the run (3.6e-11 with the reference
  # BLAS, exactly 0 with OpenBLAS 0.3.26), so only the trigger and the
  # sign of the gain are asserted; the test above asserts the maximum
  fit <- suppressWarnings(frm(ls ~ mo(income) * age,
                              data = optima_mo_data(54)))
  expect_gt(fit$opt$mo_search[["runs"]], 0)
  expect_gte(fit$opt$mo_search[["gain"]], 0)
  # and the objective's remembered best point is the fit's own
  expect_equal(unname(sdr_outer_point(fit$obj)), unname(fit$opt$par))
})

test_that("a mo() fit with random effects reports the modes of its optimum", {
  # the search evaluates the objective elsewhere; parList() reads the
  # random effects' modes from the last point evaluated, so the fit has
  # to leave the objective at its own optimum
  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300),
                   g = factor(rep(1:20, 15)))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
  fit <- frm(bf(y ~ mo(inc):z + (1 | g)) + gaussian(), data = dm)
  ref <- frmtmb:::solved_par_list(fit$obj, fit$opt$par)
  expect_equal(unname(fit$estimates$b), unname(ref$b))
})

# ---- nlminb: the point it reports ------------------------------------

optima_csmix_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}

test_that("nlminb's result is the best point it evaluated", {
  # PORT stopping on a trial it rejected ("false convergence (8)")
  # returns that trial as `par` and the best value as `objective`, and
  # frm() read its estimates from the former and restarted there. On
  # this cs() ordinal mixture (seed 3) the trial's objective is NaN.
  u <- suppressWarnings(frm(bf(y ~ x + cs(z)),
                            family = mixture(cumulative(), sratio()),
                            data = optima_csmix_data(3),
                            dry_run = "objective"))
  obj <- u$obj
  r <- suppressWarnings(frmtmb:::run_optimizer(
    "nlminb", obj$par, obj$fn, obj$gr, -Inf, Inf,
    list(eval.max = 1000, iter.max = 1000)))
  expect_true(is.finite(obj$fn(r$par)))
  expect_identical(obj$fn(r$par), r$objective)
})

# ---- cs() on an ordinal mixture --------------------------------------

test_that("a cs() ordinal mixture at its crossing wall warns, not errors", {
  # seed 3: nlminb ended on a crossed trial and the restart's first
  # gradient was NaN; seed 12: the line search converged onto two
  # thresholds that are the same double, where the objective is finite
  # and its gradient NaN. Both died with "NA/NaN gradient evaluation" on
  # 0.68.1. The optimum is on the wall where one row's category closes,
  # so the optimizer's own non-convergence is the honest report.
  for (s in c(3, 12)) {
    fit <- allow_warnings(
      frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
          data = optima_csmix_data(s)),
      c("Category specific effects", "did not report convergence",
        "Large maximum absolute gradient", "mixture(cumulative, sratio)"),
      require = "did not report convergence")
    expect_true(is.finite(fit$opt$objective))
    expect_identical(fit$obj$fn(fit$opt$par), fit$opt$objective)
  }
})

# ---- the robust log-odds in the far tails ----------------------------

test_that("a cumulative probit stays defined past |eta| = 38.2", {
  # dev/ordmix-rev-probit-sat.R: at slope 20 the 0.68.1 objective was
  # NaN, since log(pnorm()) underflowed to an infinite log-odds
  set.seed(5)
  d <- data.frame(x = rnorm(200))
  d$y <- as.integer(cut(d$x + rnorm(200), c(-Inf, -1, 0, 1, Inf)))
  fit <- frm(y ~ x, family = cumulative("probit"), data = d)
  p <- fit$opt$par
  p[names(p) == "beta"] <- 20
  expect_true(is.finite(fit$obj$fn(p)))
  expect_true(all(is.finite(fit$obj$gr(p))))
})

test_that("the cloglog and softit log-odds keep a derivative far below", {
  for (lk in c("cloglog", "softit")) {
    q <- frmtmb:::frmtmb_links[[lk]]$logit_eta
    x <- c(-1000, -745.5, -720, -50)
    tp <- RTMB::MakeTape(function(e) q(e), x)
    # eta to double precision down there
    expect_equal(tp(x), x)
    expect_identical(diag(tp$jacobian(x)), rep(1, length(x)))
  }
})

# ---- the objective's state after an escape ----------------------------

test_that("a skew normal escape reports the modes of its optimum", {
  # escape_stationary() kept the best of its runs but left the objective
  # at the last one, so parList() read the random effects' modes there:
  # seed 20 of dev/optima-rev-escape.R reported modes 0.143 away
  set.seed(20)
  g <- factor(rep(1:15, each = 10))
  u <- rnorm(15, sd = 0.8)
  d <- data.frame(g = g, x = rnorm(150))
  d$y <- 1 + 0.5 * d$x + u[g] + rnorm(150)
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = skew_normal(),
                              data = d))
  expect_false(is.null(fit$opt$stationary_escape))
  ref <- frmtmb:::solved_par_list(fit$obj, fit$opt$par)
  expect_equal(unname(fit$estimates$b), unname(ref$b))
})

# ---- a kept parameter's small loading (restores test-se-check.R) ------

test_that("a lost direction's small loading stays out of a kept row", {
  # The seed-12 mo() fixture this replaces stood on a softmax plateau
  # that lane optima removed. Built instead on a gaussian fit's own
  # objective and point, with a Hessian D S D whose S has eigenvalue
  # -0.5 along (z + 0.022 x), the seed-12 loading, and 1 elsewhere: z
  # is lost, x keeps its SE, and the null basis predictions read must be
  # zero in x's row, or every prediction through x loses its band. The
  # review's dev/optima-rev-loading2.R: without the line that zeroes a
  # kept parameter's loading that row holds 0.00189.
  set.seed(1)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- 1 + 0.5 * d$x + 0.3 * d$z + rnorm(200)
  fit <- frm(y ~ x + z, data = d)
  nm <- frmtmb:::outer_par_names(fit)
  p <- fit$opt$par
  D <- sqrt(diag(fit$obj$he(p)))
  np <- length(p)
  ix <- match("x", sub("^.*_", "", nm))
  iz <- match("z", sub("^.*_", "", nm))
  v1 <- numeric(np)
  v1[iz] <- 1
  v1[ix] <- 0.022
  v1 <- v1 / sqrt(sum(v1^2))
  Q <- qr.Q(qr(cbind(v1, diag(np))))[, seq_len(np)]
  if (sum(Q[, 1] * v1) < 0) Q[, 1] <- -Q[, 1]
  S <- Q %*% diag(c(-0.5, rep(1, np - 1))) %*% t(Q)
  r <- frmtmb:::se_tier3(fit, S * outer(D, D), matrix(0, np, np), p,
                         exact = TRUE)
  expect_true(nm[iz] %in% names(r$lost))
  expect_false(nm[ix] %in% names(r$lost))
  expect_identical(max(abs(r$null[ix, ])), 0)
})

# ---- what summary() reports ------------------------------------------

test_that("summary() reports the simplex as brms names it", {
  # the coordinates of the fit's chart have no reading of their own;
  # brms's summary lists each weight, `moincome1[1]`, under "Monotonic
  # Simplex Parameters"
  fit <- frm(ls ~ mo(income) * age, data = optima_mo_data(7))
  s <- summary(fit)
  mo <- s[["mo"]]
  expect_identical(rownames(mo), c(paste0("moincome1[", 1:3, "]"),
                                   paste0("moincome:age1[", 1:3, "]")))
  w <- mo[, "Estimate"]
  expect_equal(sum(w[1:3]), 1)
  expect_equal(sum(w[4:6]), 1)
  # a weight at a face has no delta-method error: the chart folds there
  # and the derivative is exactly 0, which printed an interval [0, 0]
  # (the review's re-check, n1); the interaction's simplex has one
  face <- w / max(w) < 1e-6
  expect_true(any(face[4:6]))
  expect_identical(attr(mo, "face"), rownames(mo)[face])
  expect_true(all(is.na(mo[face, 2:4])))
  expect_true(all(is.finite(mo[!face, "Est.Error"])))
  expect_true(all(mo[!face, 3] >= 0 & mo[!face, 4] <= 1))
  # the delta method through the fit's own covariance
  V <- vcov(fit, full = TRUE)
  j <- grep("^zeta1", rownames(V))
  z <- fit$estimates$zeta1
  J <- RTMB::MakeTape(function(x) frmtmb:::mo_simplex(x), z)$jacobian(z)
  k <- which(!face[1:3])
  expect_equal(unname(mo[k, "Est.Error"]),
               sqrt(diag(J %*% V[j, j] %*% t(J)))[k])
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("Monotonic Simplex Parameters", out, fixed = TRUE)))
  expect_true(grepl("is on the boundary of the simplex",
                    gsub("[[:space:]]+", " ", paste(out, collapse = " ")),
                    fixed = TRUE))
})

test_that("a simplex at a vertex has no delta-method interval at 1", {
  # the weight at 1 is the complement of weights at 0, so the delta
  # method degenerates there too: it gave an error near 2e-4 and printed
  # an interval of about [1, 1] (the review's final check, c1)
  fit <- frm(ls ~ mo(income) * age, data = optima_mo_data(1234))
  s <- summary(fit)
  mo <- s[["mo"]]
  w <- mo[, "Estimate"]
  tol <- frmtmb:::mo_face_tol
  vert <- which(1 - w < tol)
  expect_identical(rownames(mo)[vert], "moincome:age1[2]")
  # every weight of the simplex at the vertex, none of the interior one
  expect_identical(attr(mo, "face"), paste0("moincome:age1[", 1:3, "]"))
  expect_true(all(is.na(mo[4:6, 2:4])))
  expect_true(all(is.finite(as.matrix(mo[1:3, 2:4]))))
  out <- gsub("[[:space:]]+", " ",
              paste(utils::capture.output(print(s)), collapse = " "))
  expect_true(grepl("A weight at 0 or 1 (moincome:age1[1], moincome:age1[2]",
                    out, fixed = TRUE))
})
