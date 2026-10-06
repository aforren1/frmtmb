# An ordinal fit's thresholds and its cs() coefficients are stored under
# the names variables(fit) and brms give them, b_Intercept[k] and
# bcs_<column>[k], with brms's values: the thresholds themselves, not
# cumulative()'s first threshold and log increments. They used to be
# stored as tau_raw_k and bcs<j>_k (csfactor findings). A draws object
# stored the old way still reads, because a column under an internal
# name is on the internal scale. Against the reference build at the same
# sampler seed the stored draws, posterior_epred() and fixef() agree
# exactly (dev/sampfix-07-ordinal.R).

ord_data <- function() {
  set.seed(405L)
  n <- 300L
  do <- data.frame(x = stats::rnorm(n),
                   fc = factor(sample(c("a", "b", "c"), n, TRUE)),
                   grp = factor(sample(c("p", "q"), n, TRUE)))
  do$yo <- factor(cut(0.8 * do$x + stats::rlogis(n),
                      c(-Inf, -0.5, 0.7, Inf), labels = FALSE),
                  ordered = TRUE)
  do
}

# the ML point jittered, stored the old way (internal names and scale)
# and the new way (through the storage map)
ord_pair <- function(fit, n = 5L) {
  tpl <- fit$frame[["par_template"]]
  # a coefficient held at a constant (every ordinal family's disc) is
  # mapped, so the sampler has no column for it, as brms_par_labels()
  # has no name for it
  est <- unlist(lapply(names(tpl), function(cp) {
    v <- fit$estimates[[cp]]
    fx <- fit$frame[["betad_fixed_idx"]]
    if (identical(cp, "betad") && length(fx)) v <- v[-fx]
    v
  }))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(2)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  old <- structure(list(stanfit = NULL, draws = cbind(M, lp__ = 0),
                        fit = fit), class = "frmtmb_draws")
  new <- old
  new$draws <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  list(old = old, new = new, raw = M)
}

test_that("ordinal and cs() draws carry variables(fit)'s names", {
  do <- ord_data()
  fits <- list(
    frm(bf(yo ~ x), family = cumulative(), data = do),
    frm(bf(yo ~ x + cs(fc)), family = sratio(), data = do),
    frm(bf(yo | thres(gr = grp) ~ x), family = cumulative(), data = do))
  for (fit in fits) {
    p <- ord_pair(fit)
    expect_identical(setdiff(variables(p$new), "lp__"), variables(fit))
    expect_false(any(grepl("tau_raw|^bcs[0-9]", variables(p$new))))
    # the outer parameters that summary() and mcmc_plot() report
    expect_identical(rownames(summary(p$new)), variables(fit))
  }
})

test_that("a stored threshold is the threshold, not its internal value", {
  do <- ord_data()
  fit <- frm(bf(yo ~ x), family = cumulative(), data = do)
  p <- ord_pair(fit)
  raw <- p$raw[, c("tau_raw_1", "tau_raw_2")]
  tau <- t(apply(raw, 1L, function(r) c(r[1L], r[1L] + exp(r[2L]))))
  expect_equal(unname(p$new$draws[, c("b_Intercept[1]", "b_Intercept[2]")]),
               unname(tau))
  expect_true(all(p$new$draws[, "b_Intercept[2]"] >
                    p$new$draws[, "b_Intercept[1]"]))
})

test_that("the old storage and the new one answer every reader alike", {
  do <- ord_data()
  for (fit in list(frm(bf(yo ~ x), family = cumulative(), data = do),
                   frm(bf(yo ~ x + cs(fc)), family = sratio(), data = do))) {
    p <- ord_pair(fit)
    expect_equal(fixef(p$new), fixef(p$old))
    expect_identical(rownames(fixef(p$new)), rownames(fixef(fit)))
    expect_equal(posterior_epred(p$new), posterior_epred(p$old))
    expect_equal(posterior_linpred(p$new), posterior_linpred(p$old))
    set.seed(3)
    a <- posterior_predict(p$new)
    set.seed(3)
    b <- posterior_predict(p$old)
    expect_identical(a, b)
    h <- "Intercept[1] < Intercept[2]"
    expect_equal(hypothesis(p$new, h)$hypothesis,
                 hypothesis(p$old, h)$hypothesis)
  }
})

test_that("a threshold block with no free parameter keeps its column", {
  # sum-to-zero thresholds on a binary response: the one threshold is
  # held, so the sampler has no column for it, and brms still stores
  # b_Intercept[1] (dev/ordmix-rev-emptyblock.R). Both directions used
  # to drop it: the stored draws lacked the column, and the inverse
  # returned no column at all.
  set.seed(11)
  n <- 200L
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1L + (0.8 * d$x + stats::rlogis(n) > 0)
  expect_message(
    fit <- frm(bf(y ~ x), family = cumulative(threshold = "sum_to_zero"),
               data = d),
    "Only 2 levels")
  p <- ord_pair(fit)
  expect_identical(setdiff(variables(p$new), "lp__"), variables(fit))
  expect_identical(unname(p$new$draws[, "b_Intercept[1]"]),
                   rep(fixef(fit)["Intercept[1]", "Estimate"], nrow(p$raw)))
  back <- frmtmb.sample:::draws_to_natural(p$new$draws, fit, inverse = TRUE)
  expect_identical(back, p$old$draws)
})

test_that("frm_sample() stores the thresholds under brms's names", {
  skip_on_cran()
  skip_sampler()
  do <- ord_data()
  fit <- frm(bf(yo ~ x + cs(fc)), family = sratio(), data = do)
  ds <- allow_warnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3)),
    c("Effective Samples Size", "R-hat", "Rhat"))
  expect_identical(setdiff(variables(ds), "lp__"), variables(fit))
  # sratio's thresholds are its internal vector, so the stored draw is
  # the sampler's own value, bit for bit
  a <- rstan::extract(ds$stanfit, permuted = FALSE)
  tau <- a[, 1L, grep("^tau_raw", dimnames(a)[[3L]])]
  expect_identical(unname(ds$draws[, c("b_Intercept[1]", "b_Intercept[2]")]),
                   unname(tau))
  expect_s3_class(as_draws_rvars(ds), "draws_rvars")
  expect_identical(rownames(fixef(ds)), rownames(fixef(fit)))
})
