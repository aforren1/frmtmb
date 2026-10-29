# The me() rows of the brms log_prob identity with every prior proper
# (row 23b of test-brms-likelihood.R; dev/meme-l1.R prints the
# numbers). Function bodies only: they read helper-brms.R, which
# testthat sources before this file.
#
# frmtmb's MAP objective with a proper prior on every class is the
# posterior density in frmtmb's own coordinates theta (the latent values
# Xme themselves, log sdme, the us_chol_cor() chart of corme, log sigma,
# log sd, the group-level coefficients b). brms's log_prob with
# adjust_transform = TRUE is the same posterior in Stan's unconstrained
# coordinates u (zme, log sdme, the CPC of Lme, ...). So, with
# U(theta) = u,
#
#   -f(theta) = log_prob(U(theta)) + log |det dU/dtheta| + C
#
# where C collects normalizing constants (brms's half-Student-t
# lccdf terms). The Jacobian is taken by central differences of U,
# which is built by the translator brms_lp_check() uses, so it is
# independent of either density. The gradient identity is the chain
# rule, -grad f = J' grad log_prob + grad log |det J|, with the last
# term from the closed form `Lfun` of each case, which the run checks
# against the numeric log |det J| at every point.

me_lp_data <- function(seed = 23, n = 40) {
  set.seed(seed)
  tx <- stats::rnorm(n, 1, 0.8)
  d <- data.frame(x = tx + stats::rnorm(n, 0, 0.3), sx = 0.3,
                  g = factor(rep(seq_len(n / 5), each = 5)),
                  h = factor(rep(1:5, n / 5)))
  ng <- nlevels(d$g)
  d$xg <- stats::rnorm(ng)[d$g]
  d$zg <- 0.6 * d$xg + stats::rnorm(ng, 0, 0.6)[d$g]
  d$sg <- 0.25
  d$y <- 2 + 0.7 * tx + 0.5 * d$xg - 0.3 * d$zg +
    stats::rnorm(5, 0, 0.4)[d$h] + stats::rnorm(n, 0, 0.5)
  d
}

# one prior table, written once for each package
me_lp_priors <- function(rows) {
  fp <- NULL
  bp <- NULL
  for (r in rows) {
    a <- list(r$prior, class = r$class, coef = r$coef %||% "",
              group = r$group %||% "")
    one <- do.call(set_prior, a)
    fp <- if (is.null(fp)) one else fp + one
    bone <- do.call(brms::set_prior, a)
    bp <- if (is.null(bp)) bone else c(bp, bone)
  }
  list(frm = fp, brm = bp)
}

me_lp_cases <- function() {
  d <- me_lp_data()
  base <- list(list(prior = "normal(0, 1)", class = "b"),
               list(prior = "student_t(3, 0, 2.5)", class = "Intercept"),
               list(prior = "student_t(3, 0, 2.5)", class = "sigma"),
               list(prior = "normal(0, 2)", class = "meanme"),
               list(prior = "exponential(1)", class = "sdme"))
  n <- nrow(d)
  ng <- nlevels(d$g)
  nh <- nlevels(d$h)
  one_lsd <- function(p, cmp) p[cmp == "logsdme"]
  list(
    list(label = "one me()", data = d,
         bform = brms::bf(y ~ me(x, sx)), fform = bf(y ~ me(x, sx)),
         prior = me_lp_priors(base),
         # z = (X - meanme) / sdme
         Lfun = function(p, cmp) -n * sum(one_lsd(p, cmp))),
    list(label = "two me(gr = g), corme", data = d,
         bform = brms::bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)),
         fform = bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)),
         prior = me_lp_priors(c(base, list(
           list(prior = "normal(1, 0.5)", class = "meanme", coef = "mezg"),
           list(prior = "student_t(3, 0, 1)", class = "sdme",
                coef = "mexg"),
           list(prior = "lkj(2)", class = "corme", group = "g")))),
         # z = L^-1 D^-1 (X - meanme) per level, with diag(L) = (1,
         # 1 / sqrt(1 + t^2)); Stan's CPC of Lme is asinh(t)
         Lfun = function(p, cmp) {
           t <- p[cmp == "thetame"]
           -ng * sum(one_lsd(p, cmp)) + (ng - 1) * 0.5 * log(1 + t^2)
         }),
    list(label = "me() + (1 | h)", data = d,
         bform = brms::bf(y ~ me(x, sx) + (1 | h)),
         fform = bf(y ~ me(x, sx) + (1 | h)),
         prior = me_lp_priors(c(base, list(
           list(prior = "student_t(3, 0, 2.5)", class = "sd")))),
         # and z_1 = b / sd for the group-level block
         Lfun = function(p, cmp) {
           -n * sum(one_lsd(p, cmp)) - nh * sum(p[cmp == "theta"])
         }))
}

# The identity at `npts` points: the mode of frm()'s MAP fit and
# npts - 1 points around it (seeded). Returns the per-point table and
# the three summary numbers the test asserts on.
me_lp_points <- function(cs, npts = 3, seed = 1, h = 1e-6) {
  d <- cs$data
  fit <- frm(cs$fform, family = gaussian(), data = d, prior = cs$prior$frm,
             control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
  bp <- cs$prior$brm
  bfam <- stats::gaussian()
  code <- with_brms_me(function() {
    brms::make_stancode(cs$bform, data = d, family = bfam, prior = bp)
  })
  sdat <- with_brms_me(function() {
    brms_standata(cs$bform, data = d, family = bfam, prior = bp)
  })
  rtab <- with_brms_me(function() brms_ranef_table(cs$bform, d, bfam, bp))
  mod <- brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  env <- fit$obj$env
  p0 <- env$last.par.best
  cmp <- names(p0)
  at <- function(p) {
    f2 <- fit
    for (cp in unique(cmp)) f2$estimates[[cp]][] <- p[cmp == cp]
    f2
  }
  U <- function(p) {
    rstan::unconstrain_pars(sf, stan_pars_from_fit(at(p), sdat, code, rtab))
  }
  step <- function(p, j) {
    replace(numeric(length(p)), j, h * max(1, abs(p[j])))
  }
  if (length(U(p0)) != length(p0)) {
    stop("brms has ", length(U(p0)), " unconstrained parameters and ",
         "frmtmb ", length(p0))
  }
  set.seed(seed)
  pts <- c(list(p0), lapply(seq_len(npts - 1L), function(k) {
    p0 + stats::rnorm(length(p0), 0, 0.1)
  }))
  rows <- lapply(pts, function(p) {
    u <- U(p)
    J <- vapply(seq_along(p), function(j) {
      e <- step(p, j)
      (U(p + e) - U(p - e)) / (2 * e[j])
    }, numeric(length(p)))
    l_num <- as.numeric(determinant(J, logarithm = TRUE)$modulus)
    l_ana <- sum(cs$Lfun(p, cmp))
    ours <- -env$f(p)
    lpb <- rstan::log_prob(sf, u, adjust_transform = TRUE)
    g_ours <- -as.numeric(env$f(p, order = 1))
    g_b <- as.numeric(rstan::grad_log_prob(sf, u, adjust_transform = TRUE))
    g_l <- vapply(seq_along(p), function(j) {
      e <- step(p, j)
      (sum(cs$Lfun(p + e, cmp)) - sum(cs$Lfun(p - e, cmp))) / (2 * e[j])
    }, 0)
    g_pred <- as.numeric(crossprod(J, g_b)) + g_l
    c(ours = ours, brms = lpb, Lnum = l_num, Lana = l_ana,
      const = ours - lpb - l_num, grad_err = max(abs(g_ours - g_pred)),
      grad_scale = max(abs(g_ours)))
  })
  tab <- do.call(rbind, rows)
  list(table = tab,
       spread = diff(range(tab[, "const"])) /
         max(abs(tab[, c("ours", "brms", "Lnum")])),
       grad_rel = max(tab[, "grad_err"] / pmax(1, tab[, "grad_scale"])),
       L_err = max(abs(tab[, "Lnum"] - tab[, "Lana"]) /
                     pmax(1, abs(tab[, "Lnum"]))))
}
