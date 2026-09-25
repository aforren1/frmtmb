# brms's gr(g, by = f) and mm(g1, g2, by = ): one set of group-level
# standard deviations and correlations per level of the by-variable.
#
# frmtmb fits the term as one ordinary block per by-level (R/gr-by.R),
# so what these tests pin is the model (against lme4 and glmmTMB, which
# write one term per by-level with indicator columns, and against
# brms's density written from its Stan code), brms's names, the views
# brms keys by the grouping factor, the routing of a new level to the
# covariance of its own by-level, and each refusal. dev/grby-validate.R
# is the full validation and dev/grby-findings.md its record.
#
# Every numeric tolerance is relative to the quantity compared.

grby_data <- function(seed = 11, ng = 45, per = 8) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), each = per)))
  d$f <- factor(c("a", "b", "c")[(as.integer(d$g) - 1L) %% 3L + 1L])
  d$x <- rnorm(nrow(d))
  sd0 <- rbind(a = c(1.0, 0.2), b = c(0.3, 0.8), c = c(0.6, 0.5))
  k <- as.character(d$f[match(seq_len(ng), as.integer(d$g))])
  u <- cbind(rnorm(ng, 0, sd0[k, 1]), rnorm(ng, 0, sd0[k, 2]))
  d$y <- 1 + 0.5 * d$x + u[d$g, 1] + u[d$g, 2] * d$x +
    rnorm(nrow(d), 0, 0.5)
  for (l in c("a", "b", "c")) {
    d[[paste0("f", l)]] <- as.numeric(d$f == l)
    d[[paste0("x", l)]] <- d[[paste0("f", l)]] * d$x
  }
  d
}

grby_fit <- local({
  cache <- list()
  function(form = y ~ x + (1 + x | gr(g, by = f)), ...) {
    key <- deparse1(form)
    if (is.null(cache[[key]])) {
      cache[[key]] <<- frm(bf(form) + gaussian(), data = grby_data(), ...)
    }
    cache[[key]]
  }
})

# ------------------------------------------------------------ the model

test_that("the term is one block per by-level over that level's groups", {
  fr <- frm(bf(y ~ x + (1 + x | gr(g, by = f))) + gaussian(),
            data = grby_data(), dry_run = "frame")
  bks <- fr$re_blocks
  expect_length(bks, 3L)
  expect_identical(vapply(bks, function(bk) bk$by$level, ""),
                   c("a", "b", "c"))
  expect_identical(vapply(bks, `[[`, 0L, "n_levels"), c(15L, 15L, 15L))
  # every level of g is in exactly one block, and the block is the one
  # its rows' f names
  d <- grby_data()
  for (bk in bks) {
    rows <- d$g %in% bk$levels
    expect_true(all(d$f[rows] == bk$by$level))
  }
  expect_setequal(unlist(lapply(bks, `[[`, "levels")), levels(d$g))
})

test_that("the fit agrees with lme4 fitting one term per by-level", {
  skip_if_not_installed("lme4")
  fit <- grby_fit()
  ref <- lme4::lmer(y ~ x + (0 + fa + xa | g) + (0 + fb + xb | g) +
                      (0 + fc + xc | g), data = grby_data(), REML = FALSE,
                    control = lme4::lmerControl(
                      check.conv.singular = "ignore"))
  ll <- as.numeric(logLik(ref))
  expect_lt(abs(as.numeric(logLik(fit)) - ll), abs(ll) * 1e-8)
  vl <- lme4::VarCorr(ref)
  sd_l <- unlist(lapply(vl, function(v) attr(v, "stddev")), use.names = FALSE)
  sd_f <- VarCorr(fit)$g$sd[, "Estimate"]
  expect_lt(max(abs(sd_f / sd_l - 1)), 1e-3)
  cor_l <- vapply(vl, function(v) attr(v, "correlation")[1, 2], 1)
  cr <- VarCorr(fit)$g$cor
  cor_f <- c(cr["Intercept:fa", "Estimate", "x:fa"],
             cr["Intercept:fb", "Estimate", "x:fb"],
             cr["Intercept:fc", "Estimate", "x:fc"])
  # a correlation's own scale is its distance from +-1
  expect_lt(max(abs(cor_f - unname(cor_l)) / (1 - abs(unname(cor_l)))), 1e-3)
})

test_that("logLik is brms's marginal density at the estimates", {
  # brms: r_j = diag(sd[, Jby[j]]) L[Jby[j]] z_j (scale_r_cor_by()), so a
  # gaussian response is y ~ N(X beta, Z S Z' + sigma^2 I) with S block
  # diagonal and level j taking its by-level's covariance. The Laplace
  # approximation is exact here, so this is an identity up to rounding.
  fit <- grby_fit()
  d <- grby_data()
  vc <- VarCorr(fit)
  S <- vc$g$cov[, "Estimate", ]
  V <- diag(vc$residual__$sd[1, "Estimate"]^2, nrow(d))
  for (j in levels(d$g)) {
    i <- which(d$g == j)
    k <- paste0(c("Intercept:f", "x:f"), d$f[i[1L]])
    Zi <- cbind(1, d$x[i])
    V[i, i] <- V[i, i] + Zi %*% S[k, k] %*% t(Zi)
  }
  mu <- drop(cbind(1, d$x) %*% fixef(fit)[, "Estimate"])
  R <- chol(V)
  ll <- -sum(log(diag(R))) - 0.5 * nrow(d) * log(2 * pi) -
    0.5 * sum(backsolve(R, d$y - mu, transpose = TRUE)^2)
  expect_lt(abs(as.numeric(logLik(fit)) - ll), abs(ll) * 1e-10)
})

test_that("a structured block splits too: ar1 against glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(13)
  da <- expand.grid(t = factor(1:4), g = factor(1:30))
  da$f <- factor(ifelse(as.integer(da$g) <= 15, "a", "b"))
  ua <- t(vapply(1:30, function(j) {
    r <- if (j <= 15) 0.7 else -0.2
    (if (j <= 15) 1 else 0.5) *
      as.vector(t(chol(r^abs(outer(1:4, 1:4, "-")))) %*% rnorm(4))
  }, numeric(4)))
  da$y <- ua[cbind(as.integer(da$g), as.integer(da$t))] +
    rnorm(nrow(da), 0, 0.4)
  da$fa <- as.numeric(da$f == "a")
  da$fb <- as.numeric(da$f == "b")
  fit <- frm(bf(y ~ 1 + ar1(0 + t | gr(g, by = f))) + gaussian(), data = da)
  ref <- suppressWarnings(glmmTMB::glmmTMB(
    y ~ 1 + ar1(0 + t:fa | g) + ar1(0 + t:fb | g), data = da))
  ll <- as.numeric(logLik(ref))
  expect_lt(abs(as.numeric(logLik(fit)) - ll), abs(ll) * 1e-8)
})

test_that("mm(g1, g2, by = ) is brms's density at the estimates", {
  set.seed(3)
  n <- 300
  lv <- 1:24
  dm <- data.frame(g1 = factor(sample(lv, n, TRUE), levels = lv),
                   g2 = factor(sample(lv, n, TRUE), levels = lv))
  fl <- ifelse(lv <= 12, "a", "b")
  dm$f1 <- factor(fl[dm$g1])
  dm$f2 <- factor(fl[dm$g2])
  dm$x <- rnorm(n)
  um <- rnorm(24, 0, ifelse(lv <= 12, 1, 0.3))
  dm$y <- 0.5 + 0.3 * dm$x + 0.5 * (um[dm$g1] + um[dm$g2]) +
    rnorm(n, 0, 0.5)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))) + gaussian(),
             data = dm)
  # brms's names, read off brms 2.23.0's rename_re() for this model
  # (dev/grby-log/brmsnames.txt)
  expect_identical(variables(fit)[3:4],
                   c("sd_mmg1g2__Intercept:cbind(f1, f2)1",
                     "sd_mmg1g2__Intercept:cbind(f1, f2)2"))
  s <- VarCorr(fit)$mmg1g2$sd[, "Estimate"]
  A <- matrix(0, n, 24)
  A[cbind(seq_len(n), as.integer(dm$g1))] <- 0.5
  A[cbind(seq_len(n), as.integer(dm$g2))] <-
    A[cbind(seq_len(n), as.integer(dm$g2))] + 0.5
  V <- A %*% diag(s[ifelse(lv <= 12, 1, 2)]^2) %*% t(A) +
    diag(VarCorr(fit)$residual__$sd[1, "Estimate"]^2, n)
  mu <- drop(cbind(1, dm$x) %*% fixef(fit)[, "Estimate"])
  R <- chol(V)
  ll <- -sum(log(diag(R))) - 0.5 * n * log(2 * pi) -
    0.5 * sum(backsolve(R, dm$y - mu, transpose = TRUE)^2)
  expect_lt(abs(as.numeric(logLik(fit)) - ll), abs(ll) * 1e-10)
  expect_identical(dim(ranef(fit)[[1L]]), c(24L, 1L))
})

# ------------------------------------------------------------ the names

test_that("the parameters carry brms's names", {
  fit <- grby_fit()
  # brms 2.23.0's rename_re() for (1 + x | gr(g, by = f))
  # (dev/grby-brmsnames.R, dev/grby-log/brmsnames.txt)
  expect_identical(
    variables(fit),
    c("b_Intercept", "b_x",
      "sd_g__Intercept:fa", "sd_g__x:fa", "cor_g__Intercept:fa__x:fa",
      "sd_g__Intercept:fb", "sd_g__x:fb", "cor_g__Intercept:fb__x:fb",
      "sd_g__Intercept:fc", "sd_g__x:fc", "cor_g__Intercept:fc__x:fc",
      "sigma"))
  vc <- VarCorr(fit)
  expect_identical(names(vc), c("g", "residual__"))
  expect_identical(rownames(vc$g$sd),
                   c("Intercept:fa", "x:fa", "Intercept:fb", "x:fb",
                     "Intercept:fc", "x:fc"))
  # a pair across two by-levels is independent by construction, the 0
  # brms fills in for a correlation it does not have
  expect_identical(vc$g$cor["Intercept:fa", "Estimate", "Intercept:fb"], 0)
  rs <- summary(fit)$random$g
  expect_true(all(c("sd(Intercept:fa)", "cor(Intercept:fb,x:fb)") %in%
                    rownames(rs)))
  h <- hypothesis(fit, "sd_g__Intercept:fa - sd_g__Intercept:fb = 0",
                  class = NULL)
  sds <- vc$g$sd[, "Estimate"]
  expect_equal(h$hypothesis$Estimate,
               unname(sds["Intercept:fa"] - sds["Intercept:fb"]),
               tolerance = 1e-8)
  expect_message(ci <- confint(fit, parm = "sd_g__x:fb"), "theta_")
  expect_equal(nrow(ci), 1L)
  expect_lt(ci[1, "lwr"], ci[1, "upr"])
})

test_that("a by-split on sigma takes the dpar prefix before the by-level", {
  d <- grby_data()
  fit <- frm(bf(y ~ x, sigma ~ (1 | gr(g, by = f))) + gaussian(), data = d)
  expect_identical(grep("^sd_", variables(fit), value = TRUE),
                   paste0("sd_g__sigma_Intercept:f", c("a", "b", "c")))
})

test_that("ranef(), coef() and ngrps() have one entry over every level", {
  fit <- grby_fit()
  d <- grby_data()
  r <- ranef(fit, condVar = TRUE)
  expect_identical(names(r), "g")
  expect_identical(dimnames(r$g), list(levels(d$g), c("Intercept", "x")))
  expect_identical(dim(attr(r$g, "condSD")), c(45L, 2L))
  expect_identical(attr(r$g, "term"), "1 + x | g")
  # each level's mode is the one its own by-level's block holds
  cvec <- coef_b(fit)
  for (bk in fit$frame$re_blocks) {
    M <- t(matrix(cvec[bk$c_idx], nrow = 2L))
    expect_identical(unname(r$g[bk$levels, ]), M)
  }
  cf <- coef(fit)$g
  expect_identical(dim(cf), c(45L, 2L))
  expect_equal(cf[["x"]], unname(fixef(fit)["x", "Estimate"] + r$g[, "x"]),
               tolerance = 1e-12)
  expect_identical(ngrps(fit), list(g = 45L))
})

# ------------------------------------------------------ new data

test_that("a new level takes the covariance of its own by-level", {
  fit <- grby_fit()
  d <- grby_data()
  S <- VarCorr(fit)$g$cov[, "Estimate", ]
  nd <- data.frame(g = c("new1", "new2", "new3"), f = c("a", "b", "c"),
                   x = 0.7)
  p <- frm_linpred(fit, newdata = nd, allow_new_levels = TRUE, se.fit = TRUE)
  base <- frm_linpred(fit, newdata = nd, se.fit = TRUE, re_formula = NA)
  extra <- unname(p$se.fit^2 - base$se.fit^2)
  z <- c(1, 0.7)
  want <- vapply(c("a", "b", "c"), function(k) {
    i <- paste0(c("Intercept:f", "x:f"), k)
    drop(t(z) %*% S[i, i] %*% z)
  }, 1)
  expect_equal(extra, unname(want), tolerance = 1e-8)
  # an old level keeps its fitted effect whatever f the new row says
  old <- d[c(1L, 10L), c("g", "f", "x")]
  swapped <- old
  swapped$f <- factor(c("c", "a"), levels = levels(d$f))
  expect_identical(frm_linpred(fit, newdata = swapped),
                   frm_linpred(fit, newdata = old))
  expect_equal(as.vector(frm_linpred(fit, newdata = old)),
               unname(fitted(fit)[c(1L, 10L), "Estimate"]),
               tolerance = 1e-12)
  # and needs no by-variable at all
  expect_identical(frm_linpred(fit, newdata = old[, c("g", "x")]),
                   frm_linpred(fit, newdata = old))
})

test_that("new data that cannot be routed is refused by name", {
  fit <- grby_fit()
  expect_error(predict(fit, newdata = data.frame(g = "new", f = "a", x = 1)),
               "New levels in grouping factor `g`: new")
  expect_error(frm_linpred(fit, newdata = data.frame(g = "new", f = "zz",
                                                     x = 1),
                           allow_new_levels = TRUE),
               "fall in f = 'zz', which has no fitted covariance")
  expect_error(frm_linpred(fit, newdata = data.frame(g = "new", f = NA,
                                                     x = 1),
                           allow_new_levels = TRUE),
               "gr[(]by = f[)] is missing")
  expect_error(frm_linpred(fit, newdata = data.frame(g = "new",
                                                     f = c("a", "b"), x = 1),
                           allow_new_levels = TRUE),
               "Some levels of 'g' correspond to multiple levels of 'f'")
  # a by column newdata does not carry; its name is one nothing else
  # in scope carries, so R's lookup cannot find it elsewhere
  d <- grby_data()
  d$grby_only_col <- d$f
  fz <- frm(bf(y ~ x + (1 | gr(g, by = grby_only_col))) + gaussian(),
            data = d)
  expect_error(frm_linpred(fz, newdata = data.frame(g = "new", x = 1),
                           allow_new_levels = TRUE),
               "cannot evaluate the by-variable on newdata")
})

test_that("simulate() draws from the split blocks", {
  fit <- grby_fit()
  s <- simulate(fit, nsim = 2, seed = 1)
  expect_identical(dim(s), c(nrow(grby_data()), 2L))
  expect_false(anyNA(s))
})

test_that("frm_simulate() takes the by-level names", {
  d <- grby_data()
  np <- list(b_Intercept = 1, b_x = 0.5, sigma = 0.5,
             "sd_g__Intercept:fa" = 1, "sd_g__Intercept:fb" = 0.2,
             "sd_g__Intercept:fc" = 0.6)
  s <- frm_simulate(bf(y ~ x + (1 | gr(g, by = f))) + gaussian(), d,
                    newparams = np, nsim = 1, seed = 1)
  expect_identical(dim(s), c(nrow(d), 1L))
  expect_error(frm_simulate(bf(y ~ x + (1 | gr(g, by = f))) + gaussian(), d,
                            newparams = list(b_Intercept = 1, b_x = 0.5,
                                             sigma = 0.5,
                                             sd_g__Intercept = 1)),
               "sd_g__Intercept:fa")
})

# ------------------------------------------------------------ refusals

test_that("a level of g in two by-levels is refused with brms's message", {
  d <- grby_data()
  d$h <- factor(rep(c("u", "v"), length.out = nrow(d)))
  expect_error(frm(bf(y ~ x + (1 | gr(g, by = h))) + gaussian(), data = d),
               "Some levels of 'g' correspond to multiple levels of 'h'")
})

test_that("brms's other by refusals are brms's", {
  d <- grby_data()
  d$f2 <- d$f
  expect_error(frm(bf(y ~ x + (1 | gr(g, by = f)) +
                        (0 + x | gr(g, by = f2))) + gaussian(), data = d),
               "Each grouping factor can only be associated with one 'by'")
  set.seed(1)
  dm <- data.frame(g1 = factor(sample(1:6, 60, TRUE)),
                   g2 = factor(sample(1:6, 60, TRUE)), y = rnorm(60))
  dm$f1 <- factor(as.integer(dm$g1) <= 3)
  expect_error(frm(bf(y ~ (1 | mm(g1, g2, by = f1))) + gaussian(),
                   data = dm),
               "Grouping structure 'mm' expects 'by' to be a matrix")
})

test_that("structures that do not split are refused by name", {
  d <- grby_data()
  A <- diag(45)
  dimnames(A) <- list(levels(d$g), levels(d$g))
  expect_error(frm(bf(y ~ x + (1 | gr(g, by = f, cov = A))) + gaussian(),
                   data = d, data2 = list(A = A)),
               "levels in different by-levels stay correlated")
  expect_error(frm(bf(y ~ x + rr(0 + x | gr(g, by = f), d = 1)) +
                     gaussian(), data = d),
               "not supported for rr[(][)]")
  expect_error(frm(bf(y ~ x + (1 | gr(g, by = f, prec = A))) + gaussian(),
                   data = d, data2 = list(A = A)),
               "gr[(]g, by = , prec = [)] is not supported")
  # the term as written is what the duplicate check reads
  expect_error(frm(bf(y ~ x + (1 | gr(g, by = f)) + (1 | g)) + gaussian(),
                   data = d),
               "Duplicated group-level effects")
})

test_that("the level-to-by-level map is brms's Jby_1", {
  # brms 2.23.0's tests.standata.R data and its asserted Nby_1 and
  # Jby_1: Stan data frmtmb does not build, whose content is the
  # by-level of each grouping level, which the fitted blocks carry
  gvar <- rep(c("1A", "1B", "2A", "2B", "3A", "3B", "10", "100", "2", "3"),
              each = 10)
  g_order <- order(gvar)
  byvar <- factor(rep(c(0, 4.5, 3, 2, "x 1"), each = 20))
  set.seed(1)
  dat <- data.frame(y = rnorm(100), x = rnorm(100), g = gvar,
                    g2 = gvar[g_order], z = byvar, z2 = byvar[g_order])
  jby <- c(2, 2, 1, 1, 5, 4, 4, 5, 3, 3)
  for (form in list(y ~ x + (x | gr(g, by = z)),
                    y ~ x + (x | mm(g, g2, by = cbind(z, z2))))) {
    fr <- frm(bf(form) + gaussian(), data = dat, dry_run = "frame")
    by <- fr$re_blocks[[1L]]$by
    expect_length(fr$re_blocks, 5L)
    expect_identical(match(by$level_by, by$levels), as.integer(jby))
  }
  # and brms's names drop the level's whitespace, "x 1" is zx1
  fr <- frm(bf(y ~ x + (1 | gr(g, by = z))) + gaussian(), data = dat,
            dry_run = "frame")
  expect_identical(fr$re_blocks[[5L]]$by$name, "zx1")
})

test_that("the compatibility registry states what the fits do", {
  expect_equal(frm_compat("gr_by", "us")$status, "works")
  expect_equal(frm_compat("gr_by", "ar1")$status, "works")
  expect_equal(frm_compat("gr_by", "mm()")$status, "works")
  for (b in c("gr_cov", "gr_prec", "rr", "equalto", "importance")) {
    expect_equal(frm_compat("gr_by", b)$status, "refused")
  }
  expect_error(frm(bf(y ~ x + (1 | gr(g, by = f))) + gaussian(),
                   data = grby_data(), importance = 50L),
               "carry the same grouping levels")
})

test_that("several by-split terms on one factor are ordered as brms orders", {
  d <- grby_data()
  fit <- frm(bf(y ~ x + (1 + x || gr(g, by = f))) + gaussian(), data = d)
  # brms's get_rnames() is by-level first within the group, every term's
  # coefficients inside one by-level
  expect_identical(rownames(VarCorr(fit)$g$sd),
                   paste0(rep(c("Intercept:f", "x:f"), 3),
                          rep(c("a", "b", "c"), each = 2)))
  # summary.brmsfit() reads sd_ then cor_ rows in variables() order
  rs <- rownames(summary(grby_fit())$random$g)
  expect_identical(rs, c(paste0("sd(", c("Intercept", "x"), ":f",
                                rep(c("a", "b", "c"), each = 2), ")"),
                         paste0("cor(Intercept:f", c("a", "b", "c"), ",x:f",
                                c("a", "b", "c"), ")")))
})
