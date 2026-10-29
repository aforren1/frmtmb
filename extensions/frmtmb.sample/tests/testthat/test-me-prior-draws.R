# frm_sample() on me() terms: brms 2.23.0's default priors on meanme,
# sdme and corme, prior draws of the three classes under brms's names,
# and "only" (dev/meme-findings.md; brms's rows and names in
# dev/meme-out/brms.txt, brms-order.txt and brms-names.txt). Through
# frmtmb.sample 0.24.0 the sampler was flat on log(sdme) and on the
# chart of the correlation, where brms is flat on sdme and lkj(1) on the
# correlation matrix; the classes had no slot, so no prior could be
# written on them, "yes" drew none of them and "only" refused.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

mpd_data <- function(seed = 61, n = 60) {
  set.seed(seed)
  tx <- stats::rnorm(n, 1, 0.8)
  tz <- 0.5 * tx + stats::rnorm(n, 0, 0.7)
  tw <- 0.3 * tx + stats::rnorm(n, 0, 0.8)
  d <- data.frame(x = tx + stats::rnorm(n, 0, 0.3), sx = 0.3,
                  z = tz + stats::rnorm(n, 0, 0.3), sz = 0.3,
                  w = tw + stats::rnorm(n, 0, 0.3), sw = 0.3,
                  g = factor(rep(seq_len(n / 5), each = 5)))
  d$xg <- stats::rnorm(nlevels(d$g))[d$g]
  d$sg <- 0.25
  d$y <- 2 + 0.7 * tx - 0.4 * tz + stats::rnorm(n, 0, 0.5)
  d
}

mpd_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

mpd_ri <- function(fit, prior = NULL) {
  suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, prior))$ri
}

test_that("the sampling route lists brms's default rows for me()", {
  d <- mpd_data()
  tab <- as.data.frame(default_prior(y ~ me(x, sx) + me(z, sz) +
                                       me(xg, sg, gr = g), data = d,
                                     route = "sample"))
  me <- tab[tab$class %in% c("meanme", "sdme", "corme"), ]
  got <- stats::setNames(me$prior, paste(me$class, me$coef, me$group,
                                         sep = "|"))
  # brms 2.23.0 (dev/meme-out/brms.txt, case G): lkj(1) on corme, every
  # other row flat, sdme bounded below by 0
  want <- c("corme||" = "lkj(1)", "meanme||" = "(flat)",
            "meanme|mex|" = "(flat)", "meanme|mexg|" = "(flat)",
            "meanme|mez|" = "(flat)", "sdme||" = "(flat)",
            "sdme|mex|" = "(flat)", "sdme|mexg|" = "(flat)",
            "sdme|mez|" = "(flat)")
  expect_identical(got, want)
  expect_identical(me$lb[me$class == "sdme" & !nzchar(me$coef)], 0)
  skip_if_not_installed("brms")
  bt <- as.data.frame(local({
    attach(list(me = brms::me), name = "mpd_brms_me",
           warn.conflicts = FALSE)
    on.exit(detach("mpd_brms_me"))
    brms::default_prior(brms::bf(y ~ me(x, sx) + me(z, sz) +
                                   me(xg, sg, gr = g)), data = d)
  }))
  bt <- bt[bt$class %in% c("meanme", "sdme", "corme"), ]
  bp <- ifelse(nzchar(bt$prior), bt$prior, "(flat)")
  expect_identical(stats::setNames(bp, paste(bt$class, bt$coef, bt$group,
                                             sep = "|")), want)
})

# Flat on sdme means the sampled density on logsdme is the log map's
# Jacobian alone, and lkj(1) on corme the uniform density on the
# correlation matrices carried to the chart: -log prior + log |J| is
# the same at every point. log |J| of the correlation chart is taken by
# central differences of core's us_chol_cor(), not from the LKJ code,
# so a wrong Jacobian there cannot cancel itself here
mpd_fd_logjac <- function(M) {
  function(t) {
    f <- function(u) {
      C <- as.matrix(frmtmb::us_chol_cor(u, M))
      C[lower.tri(C)]
    }
    J <- vapply(seq_along(t), function(j) {
      h <- 1e-5 * max(1, abs(t[j]))
      e <- replace(numeric(length(t)), j, h)
      (f(t + e) - f(t - e)) / (2 * h)
    }, numeric(length(t)))
    log(abs(det(as.matrix(J))))
  }
}

test_that("the sampled default is flat on sdme and lkj(1) on corme", {
  d <- mpd_data()
  fit <- frm(y ~ me(x, sx) + me(z, sz) + me(w, sw), data = d,
             dry_run = "objective")
  ri <- mpd_ri(fit)
  comp <- vapply(ri$entries, `[[`, "", "comp")
  # meanme is held on its own scale, so brms's flat needs no entry
  expect_false("meanme" %in% comp)
  me_e <- ri$entries[comp %in% c("logsdme", "thetame")]
  expect_setequal(vapply(me_e, `[[`, "", "comp"),
                  c("logsdme", "logsdme", "logsdme", "thetame"))
  nlp <- neg_log_prior_fn(me_e)
  lj <- mpd_fd_logjac(3L)
  set.seed(3)
  pts <- replicate(6, list(ls = stats::rnorm(3), t = stats::rnorm(3)),
                   simplify = FALSE)
  parts <- vapply(pts, function(q) {
    p <- fit$frame$par_template
    p[["logsdme"]][] <- q$ls
    p[["thetame"]][] <- q$t
    c(nlp(p), sum(q$ls), lj(q$t))
  }, numeric(3))
  v <- colSums(parts)
  expect_lt(diff(range(v)) / max(abs(parts)), 1e-8)
  # and the constant is the uniform's: the volume of the 3 x 3
  # correlation matrices is pi^2 / 2
  expect_lt(abs(v[1] - log(pi^2 / 2)) / log(pi^2 / 2), 1e-7)
  # prior = "flat" is brms's (flat): still flat on sdme, and flat on the
  # Cholesky factor for corme, which is not lkj(1) at this size
  rf <- mpd_ri(fit, "flat")
  cf <- vapply(rf$entries, `[[`, "", "comp")
  th <- rf$entries[[match("thetame", cf)]]
  expect_true(isTRUE(th$dist$flat))
  expect_identical(sum(cf == "logsdme"), 3L)
})

test_that("the draws derive a grouped corme under brms's name", {
  d <- mpd_data()
  d$zg <- stats::rnorm(nlevels(d$g))[d$g]
  fit <- frm(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g) + me(x, sx) +
               me(z, sz), data = d)
  # brms 2.23.0 names them corme_g__mexg__mezg and corme__mex__mez
  # (dev/meme-out/agree-brms.txt). Through frmtmb.sample 0.24.0 the
  # grouped one was missing from the draws: the derived columns were
  # matched on "corme__" only
  nm <- frmtmb.sample:::draws_cov_cols(fit)$names
  expect_true(all(c("corme_g__mexg__mezg", "corme__mex__mez", "sdme_mexg",
                    "sdme_mex") %in% nm))
})

test_that("a prior of one's own reaches the slot it names", {
  d <- mpd_data()
  fit <- frm(y ~ me(x, sx) + me(z, sz), data = d, dry_run = "objective")
  ri <- mpd_ri(fit, set_prior("normal(0, 2)", class = "meanme") +
                 set_prior("exponential(1)", class = "sdme",
                           coef = "mez") +
                 set_prior("lkj(3)", class = "corme"))
  key <- vapply(ri$entries, function(e) paste0(e$comp, ".", e$idx[1L]), "")
  at <- function(k) ri$entries[[match(k, key)]]$dist
  expect_identical(at("meanme.1")$kind, "normal")
  expect_identical(at("meanme.2")$kind, "normal")
  # the default keeps the other sdme flat on sdme
  expect_identical(at("logsdme.1")$kind, "flat")
  expect_identical(at("logsdme.2")$kind, "exponential")
  expect_identical(at("thetame.1")$eta, 3)
})

test_that("\"yes\" draws brms's me() priors under brms's names", {
  skip_sampler()
  d <- mpd_data()
  pr <- set_prior("normal(0, 2)", class = "meanme") +
    set_prior("normal(1, 0.5)", class = "meanme", coef = "mez") +
    set_prior("exponential(1)", class = "sdme") +
    set_prior("lkj(2)", class = "corme")
  ds <- mpd_quiet(frm_sample(bf(y ~ me(x, sx) + me(z, sz) +
                                  me(xg, sg, gr = g)),
                             data = d, family = gaussian(), prior = pr,
                             sample_prior = "yes", chains = 2, iter = 1000,
                             seed = 7, refresh = 0))
  pd <- prior_draws(ds)
  # brms 2.23.0's names and order on this model and prior
  # (dev/meme-out/brms-names.txt): a coefficient's own prior splits its
  # group's vector, and brms's renaming spells that meanme_1_meanme_mex
  expect_identical(names(pd), c("Intercept", "sigma", "meanme_1_meanme_mex",
                                "meanme_1_meanme_mez", "sdme_1", "Lme_1",
                                "meanme_2", "sdme_2"))
  ks <- function(v, cdf, ...) suppressWarnings(stats::ks.test(v, cdf, ...))
  expect_gt(ks(pd$meanme_1_meanme_mex, "pnorm", 0, 2)$p.value, 0.01)
  expect_gt(ks(pd$meanme_1_meanme_mez, "pnorm", 1, 0.5)$p.value, 0.01)
  expect_gt(ks(pd$meanme_2, "pnorm", 0, 2)$p.value, 0.01)
  expect_gt(ks(pd$sdme_1, "pexp", 1)$p.value, 0.01)
  # one correlation of an LKJ(2) 2 x 2 matrix is 2 Beta(2, 2) - 1
  expect_gt(ks((pd$Lme_1 + 1) / 2, "pbeta", 2, 2)$p.value, 0.01)
  # with brms's defaults only the correlation is drawn, as in brms
  d2 <- mpd_quiet(frm_sample(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                             family = gaussian(), sample_prior = "yes",
                             chains = 1, iter = 400, seed = 7, refresh = 0))
  expect_identical(names(prior_draws(d2)), c("Intercept", "sigma", "Lme_1"))
  # draws carry brms's names for the hyperparameters
  expect_true(all(c("meanme_mex", "meanme_mez", "meanme_mexg", "sdme_mex",
                    "sdme_mez", "sdme_mexg", "corme__mex__mez") %in%
                    colnames(ds$draws)))
  ps <- as.data.frame(prior_summary(ds))
  expect_true(all(c("meanme", "sdme", "corme") %in% ps$class))
})

test_that("a bound on one coefficient is refused with brms's words", {
  d <- mpd_data()
  d$w <- stats::rnorm(nrow(d))
  fit <- frm(y ~ w + (1 + w | g) + me(x, sx) + me(z, sz), data = d,
             dry_run = "objective")
  # brms 2.23.0 stops each of these in brm() with the same sentence
  # (dev/meme-out/fix-brms-grp.txt). Through frmtmb.sample 0.24.0 they
  # sampled the bound and drew the prior without it: sdme_mex at
  # ub = 0.5 had prior draws up to 10.7 (dev/reviews/2026-10-28-meme.md,
  # finding 2)
  msg <- "Prior argument 'coef' may not be specified when using boundaries"
  # first on a model 0.24.0 could already sample, b and sd
  f0 <- frm(y ~ w + (1 + w | g), data = d, dry_run = "objective")
  for (p in list(set_prior("normal(0, 1)", class = "b", coef = "w",
                           ub = 0.1),
                 set_prior("normal(0, 1)", class = "sd", group = "g",
                           coef = "w", ub = 2))) {
    expect_error(mpd_ri(f0, p), msg, fixed = TRUE)
  }
  ri0 <- mpd_ri(f0, set_prior("normal(0, 1)", class = "sd", group = "g",
                              ub = 2))
  expect_identical(unname(ri0$upper[c("theta_1", "theta_2")]),
                   rep(log(2), 2))
  for (p in list(set_prior("normal(0, 1)", class = "b", coef = "w",
                           ub = 0.1),
                 set_prior("normal(0, 1)", class = "meanme", coef = "mex",
                           lb = 0),
                 set_prior("exponential(1)", class = "sdme", coef = "mex",
                           ub = 0.5),
                 set_prior("", class = "sdme", coef = "mex", ub = 3))) {
    expect_error(mpd_ri(fit, p), msg, fixed = TRUE)
  }
  # a class-level bound is brms's own and stays, so does a group's, and
  # frm() keeps its box on one coefficient
  ri <- mpd_ri(fit, set_prior("exponential(1)", class = "sdme", ub = 3) +
                 set_prior("normal(0, 1)", class = "sd", group = "g",
                           ub = 2) +
                 set_prior("normal(0, 1)", class = "b", lb = -1))
  expect_identical(unname(ri$upper[c("logsdme_mex", "logsdme_mez")]),
                   rep(log(3), 2))
  expect_identical(unname(ri$lower[["w"]]), -1)
  # the box binds here, which frm() reports
  fb <- allow_warnings(
    frm(y ~ w + me(x, sx), data = d,
        prior = set_prior("", class = "b", coef = "w", ub = 0.1)),
    "a bound holds it, so the fit is a constrained optimum there")
  expect_lte(fb$estimates$beta[["w"]], 0.1)
})

test_that("a class-level bound truncates the draws and the prior draws", {
  skip_sampler()
  d <- mpd_data()
  ds <- mpd_quiet(frm_sample(
    bf(y ~ me(x, sx) + me(z, sz)), data = d, family = gaussian(),
    prior = set_prior("exponential(1)", class = "sdme", ub = 0.95),
    sample_prior = "yes", chains = 2, iter = 600, seed = 2, refresh = 0))
  expect_lte(max(ds$draws[, c("sdme_mex", "sdme_mez")]), 0.95)
  ps <- prior_draws(ds)$sdme_1
  expect_lte(max(ps), 0.95)
  # truncated exponential(1) on (0, 0.95)
  cdf <- function(q) stats::pexp(q, 1) / stats::pexp(0.95, 1)
  expect_gt(suppressWarnings(stats::ks.test(ps, cdf))$p.value, 0.01)
})

# the lines of frm_sample()'s default-prior announcement
mpd_announce <- function(fit) {
  out <- character(0)
  withCallingHandlers(frmtmb.sample:::sample_resolve_priors(fit, NULL),
                      message = function(m) {
                        out <<- c(out, strsplit(conditionMessage(m),
                                                "\n")[[1L]])
                        invokeRestart("muffleMessage")
                      })
  out
}

test_that("the announcement prints delta, sdme and corme as brms does", {
  d <- mpd_data()
  fit <- frm(y ~ me(x, sx) + me(z, sz), data = d, dry_run = "objective")
  a <- mpd_announce(fit)
  expect_true(paste0("  sdme               (flat), lb = 0  ",
                     "[natural sd scale]") %in% a)
  expect_true("  corme              lkj(1)  [correlation matrix]" %in% a)
  d$yo <- cut(d$y, c(-Inf, 1, 2, 3, Inf), labels = FALSE)
  fo <- frm(yo ~ x, data = d, family = cumulative(threshold = "equidistant"),
            dry_run = "objective")
  # through frmtmb.sample 0.24.0 this line read "delta  flat()"
  expect_true("  delta              (flat), lb = 0" %in% mpd_announce(fo))
})
