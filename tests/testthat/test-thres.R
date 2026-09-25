# brms's thres() addition term on the ordinal families: thres(x = K)
# sets the number of thresholds, thres(gr = g) gives every level of g a
# threshold vector of its own. dev/thres-validate.R has the full
# comparison against brms's own R densities, MASS::polr and
# ordinal::clm; dev/thres-findings.md has the numbers.

thres_data <- function(seed = 11, n = 240) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE),
                  id = factor(sample(1:20, n, TRUE)))
  u <- rlogis(n, 0.7 * d$x + rnorm(20, 0, 0.4)[d$id])
  tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6),
              c = c(-1.5, -0.4, 0.5, 1.6))
  d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]),
                1L)
  d
}

# The category probability of each ordinal family, written from brms's
# Stan functions (cumulative_logit_lpmf and friends) and not from
# frmtmb's code: one row, its own thresholds `th`, distribution `Fd`.
thres_ref_prob <- function(fam, y, eta, th, Fd) {
  K <- length(th) + 1L
  p <- switch(fam,
    cumulative = diff(c(0, Fd(th - eta), 1)),
    sratio = {
      h <- Fd(th - eta)
      c(h, 1) * c(1, cumprod(1 - h))
    },
    cratio = {
      h <- 1 - Fd(eta - th)
      c(h, 1) * c(1, cumprod(1 - h))
    },
    acat = {
      e <- c(0, cumsum(eta - th))
      exp(e - max(e)) / sum(exp(e - max(e)))
    })
  stopifnot(length(p) == K)
  p[y]
}

thres_pt <- list(a = c(-0.8, 0.1, 1.1), b = c(-0.3, 0.9),
                 c = c(-1.2, -0.3, 0.4, 1.9))
thres_raw <- function(th, ordered) {
  unlist(lapply(th, function(t) if (ordered) c(t[1], log(diff(t))) else t))
}

test_that("grouped thresholds reproduce brms's density at a shared point", {
  d <- thres_data()
  Fs <- list(logit = stats::plogis, probit = stats::pnorm,
             cauchit = stats::pcauchy)
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    links <- switch(fam, cumulative = names(Fs), acat = "logit",
                    c("logit", "probit"))
    for (lk in links) {
      ff <- get(fam, envir = asNamespace("frmtmb"))
      fit <- frm(y | thres(gr = g) ~ x, data = d, family = ff(link = lk),
                 dry_run = "objective")$obj
      par <- fit$par
      par[names(par) == "beta"] <- 0.55
      par[names(par) == "tau_raw"] <-
        thres_raw(thres_pt, fam %in% c("cumulative", "sratio"))
      rows <- vapply(seq_len(nrow(d)), function(i) {
        log(thres_ref_prob(fam, d$y[i], 0.55 * d$x[i],
                           thres_pt[[d$g[i]]], Fs[[lk]]))
      }, 0)
      # a sum of 240 terms, each exact to a few ulps
      expect_lt(abs(-fit$fn(par) - sum(rows)),
                64 * .Machine$double.eps * sum(abs(rows)),
                label = paste(fam, lk))
    }
  }
})

test_that("grouped thresholds with group slopes are polr fitted per group", {
  skip_if_not_installed("MASS")
  d <- thres_data(seed = 12, n = 600)
  fit <- frm(y | thres(gr = g) ~ g:x, data = d, family = cumulative())
  refs <- lapply(c("a", "b", "c"), function(gg) {
    MASS::polr(factor(y) ~ x, data = d[d$g == gg, ], Hess = TRUE)
  })
  ll_ref <- sum(vapply(refs, function(r) as.numeric(logLik(r)), 0))
  ll <- as.numeric(logLik(fit))
  expect_lt(abs(ll - ll_ref) / abs(ll_ref), sqrt(.Machine$double.eps))
  fe <- fixef(fit)
  th <- fe[grepl("^Intercept", rownames(fe)), , drop = FALSE]
  th_ref <- unlist(lapply(refs, `[[`, "zeta"))
  se_ref <- unlist(lapply(refs, function(r) {
    sqrt(diag(solve(r$Hessian)))[names(r$zeta)]
  }))
  # both optimizers stop within a small fraction of a standard error
  expect_lt(max(abs(th[, "Estimate"] - th_ref) / se_ref), 1e-3)
  expect_identical(rownames(th),
                   c(paste0("Intercept[a,", 1:3, "]"),
                     paste0("Intercept[b,", 1:2, "]"),
                     paste0("Intercept[c,", 1:4, "]")))
})

test_that("a shared slope with grouped thresholds is clm(nominal = ~ g)", {
  skip_if_not_installed("ordinal")
  set.seed(13)
  n <- 600
  d <- data.frame(x = rnorm(n), g = sample(c("a", "b"), n, TRUE))
  u <- rlogis(n, 0.8 * d$x)
  d$y <- ifelse(d$g == "a", 1 + (u > -1) + (u > 0.3) + (u > 1.4),
                1 + (u > -0.2) + (u > 0.5) + (u > 2))
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
  ref <- ordinal::clm(factor(y) ~ x, nominal = ~ g, data = d)
  ll_ref <- as.numeric(logLik(ref))
  expect_lt(abs(as.numeric(logLik(fit)) - ll_ref) / abs(ll_ref),
            sqrt(.Machine$double.eps))
  se_x <- sqrt(diag(vcov(ref)))[["x"]]
  expect_lt(abs(fixef(fit)["x", "Estimate"] - coef(ref)[["x"]]) / se_x,
            1e-3)
})

test_that("thres(x = K) keeps unobserved top categories", {
  d <- thres_data(seed = 14)
  d$y <- pmin(d$y, 3)
  pr <- set_prior("normal(0, 3)", class = "Intercept")
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    ff <- get(fam, envir = asNamespace("frmtmb"))
    fit <- frm(y | thres(5) ~ x, data = d, family = ff(), prior = pr)
    expect_length(fit$estimates$tau_raw, 5L)
    expect_identical(dim(fitted(fit))[3L], 6L)
    th <- c(-1, -0.2, 0.7, 1.5, 2.4)
    ll <- sum(frmtmb:::row_lpdf(
      fit$spec$responses$y$family, d$y, d$y, list(mu = 0.55 * d$x), list(),
      list(tau_raw = thres_raw(list(th),
                               fam %in% c("cumulative", "sratio")))))
    rows <- vapply(seq_len(nrow(d)), function(i) {
      log(thres_ref_prob(fam, d$y[i], 0.55 * d$x[i], th, stats::plogis))
    }, 0)
    expect_lt(abs(ll - sum(rows)), 64 * .Machine$double.eps * sum(abs(rows)),
              label = fam)
  }
  # without a prior the thresholds above category 3 are not placed
  expect_warning(frm(y | thres(5) ~ x, data = d, family = cumulative()),
                 "Intercept\\[3\\], Intercept\\[4\\], Intercept\\[5\\]")
  expect_no_warning(frm(y | thres(5) ~ x, data = d, family = cumulative(),
                        prior = pr))
  # cs() takes one coefficient per threshold, so four here
  fcs <- suppressWarnings(frm(y | thres(4) ~ cs(x), data = d,
                              family = sratio(), prior = pr))
  expect_length(fcs$estimates[[grep("^bcs", names(fcs$estimates),
                                     value = TRUE)]], 4L)
})

test_that("a per-level count in thres(x, gr) is brms's", {
  d <- thres_data()
  d$nth <- c(a = 3, b = 3, c = 4)[d$g]
  fit <- frm(y | thres(nth, g) ~ x, data = d, family = acat(),
             prior = set_prior("normal(0, 3)", class = "Intercept"))
  th <- fit$spec$responses$y$family[["thres"]]
  expect_identical(th$nthres, c(3L, 3L, 4L))
  expect_true("b_Intercept[b,3]" %in% variables(fit))
})

test_that("threshold priors are per group, uncentered, and take group =", {
  d <- thres_data()
  f0 <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative(),
            dry_run = "objective")$obj
  f1 <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative(),
            prior = set_prior("normal(0, 2)", class = "Intercept"),
            dry_run = "objective")$obj
  f2 <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative(),
            prior = set_prior("normal(0, 2)", class = "Intercept",
                              group = "b"),
            dry_run = "objective")$obj
  par <- f0$par
  par[names(par) == "beta"] <- 0.55
  par[names(par) == "tau_raw"] <- thres_raw(thres_pt, TRUE)
  # each level is an ordered vector of its own: the density is on its
  # thresholds, plus the log-Jacobian of its own increments, with no
  # centering offset
  lp <- function(t) sum(dnorm(t, 0, 2, log = TRUE)) + sum(log(diff(t)))
  want1 <- sum(vapply(thres_pt, lp, 0))
  got1 <- f0$fn(par) - f1$fn(par)
  expect_lt(abs(got1 - want1), 64 * .Machine$double.eps * abs(f0$fn(par)))
  got2 <- f0$fn(par) - f2$fn(par)
  expect_lt(abs(got2 - lp(thres_pt$b)),
            64 * .Machine$double.eps * abs(f0$fn(par)))
  expect_error(frm(y | thres(gr = g) ~ x, data = d, family = cumulative(),
                   prior = set_prior("normal(0, 2)", class = "Intercept",
                                     group = "q")),
               "whose levels are")
})

test_that("group = on an ungrouped threshold prior is refused", {
  # this used to be applied to the whole threshold vector, silently
  d <- thres_data()
  expect_error(frm(y ~ x, data = d, family = cumulative(),
                   prior = set_prior("normal(0, 2)", class = "Intercept",
                                     group = "a")),
               "one threshold vector. Drop group")
})

test_that("fitted, predict and simulate read each row's own thresholds", {
  d <- thres_data()
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = sratio())
  P <- frm_linpred(fit, type = "response")
  expect_identical(dim(P), c(nrow(d), 5L))
  top <- c(a = 4, b = 3, c = 5)[d$g]
  past <- outer(top, seq_len(5), `<`)
  expect_true(all(P[past] == 0))
  expect_lt(max(abs(rowSums(P) - 1)), 64 * .Machine$double.eps)
  nd <- d[c(3, 1, 7), ]
  Pn <- frm_linpred(fit, newdata = nd, type = "response")
  expect_identical(unname(Pn), unname(P[c(3, 1, 7), ]))
  fe <- fitted(fit, newdata = nd)
  expect_identical(dim(fe), c(3L, 4L, 5L))
  sims <- simulate(fit, nsim = 20, seed = 1)
  ymax <- apply(as.matrix(as.data.frame(lapply(sims, as.integer))), 1, max)
  expect_true(all(ymax <= top))
  sn <- simulate(fit, nsim = 20, seed = 1, newdata = nd)
  expect_true(all(apply(as.matrix(as.data.frame(lapply(sn, as.integer))),
                        1, max) <= top[c(3, 1, 7)]))
  pr <- predict(fit, newdata = nd, ndraws = 50)
  expect_true(all(pr[, 5] == 0 | top[c(3, 1, 7)] == 5))
  expect_error(fitted(fit, newdata = data.frame(x = 0)),
               "newdata needs that variable")
  expect_error(fitted(fit, newdata = data.frame(x = 0, g = "q")),
               "level\\(s\\) 'q'")
})

test_that("an ordered-factor response keeps its labels per group", {
  d <- thres_data()
  lv <- c("lo", "ml", "mid", "mh", "hi")
  d$f <- factor(lv[d$y], levels = c(lv, "top"), ordered = TRUE)
  fit <- frm(f | thres(gr = g) ~ x, data = d, family = cumulative())
  fi <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
  expect_identical(logLik(fit), logLik(fi))
  # "top" occurs nowhere, so it is no category, as in brms
  expect_identical(colnames(frm_linpred(fit, type = "response")), lv)
  s <- simulate(fit, nsim = 1, seed = 2)[[1L]]
  expect_identical(levels(s), lv)
  expect_error(frm(f | thres(5) ~ x, data = d, family = cumulative()),
               "ordered factor with 5 levels in the data")
})

test_that("names, hypothesis and the brms-shaped tables carry the group", {
  d <- thres_data()
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
  v <- variables(fit)
  expect_true(all(c("b_Intercept[a,3]", "b_Intercept[b,2]",
                    "b_Intercept[c,4]") %in% v))
  expect_false(any(c("b_Intercept[b,3]", "b_Intercept[4]") %in% v))
  h <- hypothesis(fit, "Intercept[c,1] < Intercept[a,1]")
  fe <- fixef(fit)
  expect_equal(h$hypothesis$Estimate,
               fe["Intercept[c,1]", "Estimate"] -
                 fe["Intercept[a,1]", "Estimate"])
  expect_identical(rownames(vcov(fit)), rownames(fe))
})

test_that("random effects, mo(), weights and emmeans compose", {
  d <- thres_data()
  fre <- frm(y | thres(gr = g) ~ x + (1 | id), data = d,
             family = cumulative())
  expect_identical(fre$opt$convergence, 0L)
  expect_identical(dim(fitted(fre))[3L], 5L)
  d$z <- rep_len(0:3, nrow(d))
  fmo <- frm(y | thres(gr = g) ~ x + mo(z), data = d, family = cratio())
  expect_true("moz" %in% rownames(fixef(fmo)))
  d$w <- 2
  fw <- frm(y | thres(gr = g) + weights(w) ~ x, data = d, family = acat())
  fd <- frm(y | thres(gr = g) ~ x, data = rbind(d, d), family = acat())
  expect_lt(abs(as.numeric(logLik(fw)) - as.numeric(logLik(fd))) /
              abs(as.numeric(logLik(fd))), sqrt(.Machine$double.eps))
  skip_if_not_installed("emmeans")
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
  em <- as.data.frame(emmeans::emmeans(fit, ~ x))
  # the latent predictor at the mean of x: no threshold enters it
  expect_equal(em$emmean, mean(d$x) * fixef(fit)["x", "Estimate"])
})

test_that("each thres() refusal names the reason", {
  d <- thres_data()
  expect_error(frm(y | thres(gr = g) ~ x, data = d, family = gaussian()),
               "not a valid addition term for family 'gaussian'")
  expect_error(frm(y | thres(gr = x) ~ x, data = d, family = cumulative()),
               "needs to be factor-like")
  d$nth <- rep_len(3:4, nrow(d))
  expect_error(frm(y | thres(nth) ~ x, data = d, family = cumulative()),
               "needs to be a single value")
  expect_error(frm(y | thres(nth, g) ~ x, data = d, family = cumulative()),
               "should be unique for each group")
  expect_error(frm(y | thres(2.5) ~ x, data = d, family = cumulative()),
               "must be a positive integer")
  expect_error(frm(y | thres(2, gr = g) ~ x, data = d,
                   family = cumulative()),
               "reaches category 4 in level 'a'")
  d1 <- d
  d1$y[d1$g == "b"] <- 1L
  expect_error(frm(y | thres(gr = g) ~ x, data = d1, family = cumulative()),
               "level\\(s\\) 'b'")
  expect_error(frm(y | thres(gr = g) ~ cs(x), data = d, family = sratio()),
               "category specific effects in models with multiple")
  expect_error(frm(y | thres(4) + thres(gr = g) ~ x, data = d,
                   family = cumulative()),
               "Duplicated addition term")
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
  expect_error(residuals(fit, type = "osa"), "not available with grouped")
})
