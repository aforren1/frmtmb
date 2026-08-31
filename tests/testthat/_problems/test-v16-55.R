# Extracted from test-v16.R:55

# test -------------------------------------------------------------------------
set.seed(61)
n_site <- 60
S <- 6
L <- matrix(0, S, 2)
L[, 1] <- c(1.2, 0.8, 0.5, -0.4, 0.9, 0.2)
L[2:S, 2] <- c(0.7, -0.6, 0.4, 0.3, -0.5)
f <- matrix(rnorm(n_site * 2), n_site, 2)
u <- f %*% t(L)
dd <- data.frame(
    y = rpois(n_site * S, exp(0.5 + as.vector(t(u)))),
    spp = factor(rep(seq_len(S), n_site)),
    site = factor(rep(seq_len(n_site), each = S))
  )
fit <- suppressWarnings(
    frm(bf(y ~ 1 + rr(spp + 0 | site, d = 2)) + poisson(), data = dd)
  )
if (requireNamespace("glmmTMB", quietly = TRUE)) {
    gt <- try(suppressWarnings(
      glmmTMB::glmmTMB(y ~ 1 + rr(spp + 0 | site, d = 2), data = dd,
                       family = poisson())
    ), silent = TRUE)
    if (!inherits(gt, "try-error") &&
        is.finite(as.numeric(logLik(gt)))) {
      expect_loglik_equal(fit, gt, tol = 1e-4)
    }
  }
V <- VarCorr(fit)[[1]]
expect_equal(sum(eigen(V, only.values = TRUE)$values > 1e-8), 2L)
d2 <- dd[dd$spp %in% c("1", "2", "3"), , drop = FALSE]
d2$spp <- droplevels(d2$spp)
fr <- suppressWarnings(
    frm(bf(y ~ 1 + rr(spp + 0 | site, d = 3)) + poisson(), data = d2)
  )
fu <- suppressWarnings(
    frm(bf(y ~ 1 + us(spp + 0 | site)) + poisson(), data = d2)
  )
expect_loglik_equal(fr, fu, tol = 1e-4)
expect_length(fitted(fit), nrow(dd))
expect_equal(unname(predict(fit, newdata = dd[1:6, ],
                              type = "response")),
               unname(fitted(fit)[1:6]), tolerance = 1e-10)
r <- ranef(fit)
expect_equal(dim(r[[1]]), c(60L, 6L))
expect_equal(nrow(as.data.frame(VarCorr(fit))), 6L + 15L)
s <- simulate(fit, nsim = 2, re.form = NA)
expect_equal(nrow(s), nrow(dd))
expect_error(predict(fit, se.fit = TRUE), "rr")
