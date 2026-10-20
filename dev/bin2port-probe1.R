# Probe 1: the categorical-like adapter against brms 2.23.0's own
# log_lik_<family>() on brms's preps (the blocks of tests.log_lik.R:481
# and tests.posterior_predict.R:330), and whether brms reads prep$lncor.
# Log: dev/bin2port-out/probe1.txt. Seed 495.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/brmsport-bin2-adapters.R")
set.seed(495)
ns <- 50
nobs <- 8
ncat <- 3
cmp <- function(label, a, b) {
  d <- abs(a - b)
  cat(sprintf("%-34s dims %s / %s  max abs diff %.3g  max rel %.3g  NA %d/%d\n",
              label, paste(dim(as.matrix(a)), collapse = "x"),
              paste(dim(as.matrix(b)), collapse = "x"), max(d, na.rm = TRUE),
              max(d / abs(b), na.rm = TRUE), sum(is.na(a)), sum(is.na(b))))
}
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(
  mu1 = array(rnorm(ns * nobs), dim = c(ns, nobs)),
  mu2 = array(rnorm(ns * nobs), dim = c(ns, nobs))
)
prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
prep$refcat <- 1
bprep <- prep
bprep$family <- brms::categorical()
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "categorical", prep = prep)
b <- sapply(1:nobs, brms:::log_lik_categorical, prep = bprep)
cmp("categorical log_lik", a, b)

prep$data$Y <- matrix(sample(1:20, nobs * ncat, TRUE), nrow = nobs)
prep$data$trials <- sample(1:20, nobs)
bprep <- prep
bprep$family <- brms::multinomial()
b <- sapply(1:nobs, brms:::log_lik_multinomial, prep = bprep)
a <- tryCatch(sapply(1:nobs, brms_adapt_log_lik, fun = "multinomial",
                     prep = prep), error = conditionMessage)
cat("multinomial log_lik (bare constructor):", if (is.character(a)) a else
      "ran", "\n")
# the same density with K given, which is what catmulti's bare
# multinomial() would resolve to
fm <- frmtmb::multinomial(K = ncat)
a <- sapply(1:nobs, function(i) {
  y <- matrix(prep$data$Y[i, ], ns, ncat, byrow = TRUE)
  as.numeric(frmtmb::row_lpdf(fm, y, y, list(mu2 = prep$dpars$mu1[, i],
                                             mu3 = prep$dpars$mu2[, i]),
                              list(), NULL))
})
cmp("multinomial(K = 3) row_lpdf", a, b)

prep$data$Y <- prep$data$Y / rowSums(prep$data$Y)
prep$dpars$phi <- rexp(ns, 10)
bprep <- prep
bprep$family <- brms::dirichlet()
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "dirichlet", prep = prep)
b <- sapply(1:nobs, brms:::log_lik_dirichlet, prep = bprep)
cmp("dirichlet log_lik", a, b)
# refcat 2
prep2 <- prep
prep2$refcat <- 2
bprep2 <- bprep
bprep2$refcat <- 2
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "dirichlet", prep = prep2)
b <- sapply(1:nobs, brms:::log_lik_dirichlet, prep = bprep2)
cmp("dirichlet log_lik, refcat 2", a, b)

prep$dpars$mu1 <- rexp(ns, 10)
prep$dpars$mu2 <- rexp(ns, 10)
prep$dpars$mu3 <- rexp(ns, 10)
bprep <- prep
bprep$family <- brms::brmsfamily("dirichlet2")
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "dirichlet2", prep = prep)
b <- sapply(1:nobs, brms:::log_lik_dirichlet2, prep = bprep)
cmp("dirichlet2 log_lik", a, b)

prep$dpars <- list(mu2 = rnorm(ns), mu3 = rnorm(ns),
                   sigma2 = rexp(ns, 10), sigma3 = rexp(ns, 10))
prep$lncor <- rbeta(ns, 2, 1)
bprep <- prep
bprep$family <- brms::brmsfamily("logistic_normal")
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "logistic_normal", prep = prep)
b <- sapply(1:nobs, brms:::log_lik_logistic_normal, prep = bprep)
cmp("logistic_normal log_lik", a, b)
# brms does not read prep$lncor: its value is the same without it
bprep0 <- bprep
bprep0$lncor <- NULL
b0 <- sapply(1:nobs, brms:::log_lik_logistic_normal, prep = bprep0)
cat("brms with and without prep$lncor identical:", identical(b, b0), "\n")
# and does read prep$dpars$lncor
bprep1 <- bprep0
bprep1$dpars$lncor <- prep$lncor
b1 <- sapply(1:nobs, brms:::log_lik_logistic_normal, prep = bprep1)
cat("brms with prep$dpars$lncor differs: max abs change",
    format(max(abs(b1 - b0)), digits = 3), "\n")
prep1 <- prep
prep1$lncor <- NULL
prep1$dpars$lncor <- prep$lncor
a1 <- sapply(1:nobs, brms_adapt_log_lik, fun = "logistic_normal",
             prep = prep1)
cmp("logistic_normal, prep$dpars$lncor", a1, b1)
# three log ratios (K = 4), three correlations per draw
prep4 <- structure(list(ndraws = ns, nobs = 2), class = "brmsprep")
prep4$dpars <- list(mu2 = rnorm(ns), mu3 = rnorm(ns), mu4 = rnorm(ns),
                    sigma2 = rexp(ns, 3), sigma3 = rexp(ns, 3),
                    sigma4 = rexp(ns, 3))
R0 <- matrix(c(1, .3, -.2, .3, 1, .4, -.2, .4, 1), 3)
prep4$dpars$lncor <- matrix(rep(c(R0[2, 1], R0[3, 1], R0[3, 2]), each = ns),
                            ns) + matrix(runif(3 * ns, -.05, .05), ns)
y4 <- matrix(rexp(8), 2)
prep4$data <- list(Y = y4 / rowSums(y4), ncat = 4)
prep4$refcat <- 1
bprep4 <- prep4
bprep4$family <- brms::brmsfamily("logistic_normal")
a <- sapply(1:2, brms_adapt_log_lik, fun = "logistic_normal", prep = prep4)
b <- sapply(1:2, brms:::log_lik_logistic_normal, prep = bprep4)
cmp("logistic_normal K = 4, 3 cors", a, b)

# four log ratios (K = 5), six DISTINCT correlations per draw (review
# 2026-10-14, finding 6): for D <= 3 the row-major and column-major
# orders of the lower triangle coincide, so only D >= 4 tells brms's
# fill order (get_cov_matrix(): for i in 2..D, for j in 1..i-1) from
# the other. The control hands the adapter the same correlations in
# column-major order; it must differ from brms.
R5 <- matrix(c(1, .30, -.20, .10, .30, 1, .25, -.15, -.20, .25, 1, .05,
               .10, -.15, .05, 1), 4)
row_major <- unlist(lapply(2:4, function(i) R5[i, seq_len(i - 1)]))
col_major <- unlist(lapply(1:3, function(j) R5[(j + 1):4, j]))
cat("row-major", format(row_major), "\ncol-major", format(col_major), "\n")
prep5 <- structure(list(ndraws = ns, nobs = 2), class = "brmsprep")
prep5$dpars <- list(mu2 = rnorm(ns), mu3 = rnorm(ns), mu4 = rnorm(ns),
                    mu5 = rnorm(ns), sigma2 = rexp(ns, 3),
                    sigma3 = rexp(ns, 3), sigma4 = rexp(ns, 3),
                    sigma5 = rexp(ns, 3))
y5 <- matrix(rexp(10), 2)
prep5$data <- list(Y = y5 / rowSums(y5), ncat = 5)
prep5$refcat <- 1
prep5$dpars$lncor <- matrix(rep(row_major, each = ns), ns)
bprep5 <- prep5
bprep5$family <- brms::brmsfamily("logistic_normal")
b <- sapply(1:2, brms:::log_lik_logistic_normal, prep = bprep5)
a <- sapply(1:2, brms_adapt_log_lik, fun = "logistic_normal", prep = prep5)
cmp("logistic_normal K = 5, 6 cors", a, b)
prep5c <- prep5
prep5c$dpars$lncor <- matrix(rep(col_major, each = ns), ns)
a <- sapply(1:2, brms_adapt_log_lik, fun = "logistic_normal", prep = prep5c)
cmp("control: column-major fill (must differ)", a, b)

## posterior_predict, the pp block's prep (tests.posterior_predict.R:330)
set.seed(1234)
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(
  mu1 = array(rnorm(ns * nobs, 0, 0.1), dim = c(ns, nobs)),
  mu2 = array(rnorm(ns * nobs, 0, 0.1), dim = c(ns, nobs))
)
prep$data <- list(Y = rep(1:ncat, 2), ncat = ncat)
prep$refcat <- 1
p <- sapply(1:nobs, brms_adapt_posterior_predict, fun = "categorical",
            prep = prep)
cat("pp categorical dim", dim(p), "range", range(p), "\n")
prep$dpars$phi <- rexp(ns, 1)
p <- brms_adapt_posterior_predict("dirichlet", 3, prep)
cat("pp dirichlet dim", dim(p), "max |rowSums - 1|",
    format(max(abs(rowSums(p) - 1)), digits = 3), "\n")
prep$dpars$mu1 <- rexp(ns, 10)
prep$dpars$mu2 <- rexp(ns, 10)
prep$dpars$mu3 <- rexp(ns, 10)
p <- brms_adapt_posterior_predict("dirichlet2", 3, prep)
cat("pp dirichlet2 dim", dim(p), "max |rowSums - 1|",
    format(max(abs(rowSums(p) - 1)), digits = 3), "\n")
prep$dpars <- list(mu2 = rnorm(ns), mu3 = rnorm(ns),
                   sigma2 = rexp(ns, 10), sigma3 = rexp(ns, 10))
prep$lncor <- rbeta(ns, 2, 1)
p <- brms_adapt_posterior_predict("logistic_normal", 3, prep)
cat("pp logistic_normal dim", dim(p), "max |rowSums - 1|",
    format(max(abs(rowSums(p) - 1)), digits = 3), "\n")
# moments of 20000 logistic-normal draws against brms's own simulator
# at one parameter point, a check that the chart and the refcat map
# draw the same distribution (not part of the ledger)
big <- structure(list(ndraws = 20000, nobs = 1), class = "brmsprep")
big$dpars <- list(mu2 = rep(0.4, 20000), mu3 = rep(-0.3, 20000),
                  sigma2 = rep(0.8, 20000), sigma3 = rep(0.5, 20000),
                  lncor = rep(0.6, 20000))
big$data <- list(Y = 1, ncat = 3)
big$refcat <- 1
pf <- brms_adapt_posterior_predict("logistic_normal", 1, big)
bbig <- big
bbig$family <- brms::brmsfamily("logistic_normal")
pb <- brms:::posterior_predict_logistic_normal(1, bbig)
cat("LN 20000 draws: frmtmb colMeans", format(colMeans(pf), digits = 4),
    " brms", format(colMeans(pb), digits = 4), "\n")
cat("LN cor of log ratios: frmtmb",
    format(cor(log(pf[, 2] / pf[, 1]), log(pf[, 3] / pf[, 1])), digits = 3),
    " brms", format(cor(log(pb[, 2] / pb[, 1]), log(pb[, 3] / pb[, 1])),
                    digits = 3), "\n")
