# Reviewer: the cov = FALSE ARMA fill of a missing newdata response
# (claim 3) against brms 2.23.0's own .predictor_arma(), run at the
# frmtmb fit's parameters with a minimal brmsprep. Data seed 31 (the
# lane's construction); RNG seeds given per block.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
pred_arma <- get(".predictor_arma", asNamespace("brms"))
set.seed(31)
G <- 30; Tn <- 8
d <- expand.grid(t = 1:Tn, g = factor(1:G))
d$x <- rnorm(nrow(d))
e <- as.vector(apply(matrix(rnorm(G * Tn), Tn, G), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$y <- 1 + 0.5 * d$x + e

fits <- list(
  ar1 = frm(bf(y ~ x + ar(t, g, p = 1)), data = d),
  ar2 = frm(bf(y ~ x + ar(t, g, p = 2)), data = d),
  arma11 = frm(bf(y ~ x + arma(t, g, p = 1, q = 1)), data = d),
  ma1 = frm(bf(y ~ x + ma(t, g, q = 1)), data = d),
  st = frm(bf(y ~ x + ar(t, g, p = 1)), data = d, family = student()))
coefs <- function(fit) {
  ac <- fit$frame$autocor[[1]]
  th <- fit$estimates$thetaac[ac$theta_idx]
  cf <- frmtmb:::autocor_cond_coefs(th, ac)
  b <- fixef(fit)[, "Estimate"]
  list(ar = cf$ar, ma = cf$ma, b0 = b[["Intercept"]], bx = b[["x"]],
       sigma = unname(sigma(fit)[1]),
       nu = if ("nu" %in% names(fit$spec$responses[[1]]$dpars)) {
         as.numeric(frm_linpred(fit, dpar = "nu", type = "response"))[1]
       })
}
jlag <- function(gr, max_lag) {
  N <- length(gr); out <- rep(0, N)
  for (n in seq_len(N)[-N]) {
    ind <- n:max(1, n + 1 - max_lag)
    out[n] <- sum(gr[ind] %in% gr[n + 1])
  }
  out
}
brms_fill <- function(cf, nd, S, fam = "gaussian", nu = NULL) {
  mu <- cf$b0 + cf$bx * nd$x
  eta <- matrix(mu, S, length(mu), byrow = TRUE)
  prep <- structure(list(family = list(fun = fam),
                         dpars = list(sigma = cf$sigma, nu = nu),
                         ndraws = S, data = list()), class = "brmsprep")
  ml <- max(length(cf$ar), length(cf$ma), 1)
  pred_arma(eta, ar = if (length(cf$ar)) matrix(cf$ar, S, length(cf$ar),
                                                byrow = TRUE),
            ma = if (length(cf$ma)) matrix(cf$ma, S, length(cf$ma),
                                           byrow = TRUE),
            Y = nd$y, J_lag = jlag(as.character(nd$g), ml), fprep = prep)
}
frm_fill <- function(fit, nd, cf) {
  rspec <- fit$spec$responses[[1]]
  mu0 <- cf$b0 + cf$bx * nd$x
  dfn <- function(f) {
    out <- list(mu = mu0, sigma = rep(cf$sigma, nrow(nd)))
    if (!is.null(cf$nu)) out$nu <- rep(cf$nu, nrow(nd))
    out
  }
  frmtmb:::arma_cond_fill_dpars(fit, rspec, nd, dfn)$mu
}
cat("\n== 1. exact: one group, same seed, S = 1 ==\n")
cases <- list(
  run4 = function(nd) { nd$y[nd$t >= 5] <- NA; nd },
  first_and_3 = function(nd) { nd$y[nd$t %in% c(1, 3)] <- NA; nd },
  all = function(nd) { nd$y <- NA_real_; nd },
  alt = function(nd) { nd$y[nd$t %% 2 == 0] <- NA; nd })
for (fn in names(fits)) {
  fit <- fits[[fn]]; cf <- coefs(fit)
  fam <- if (fn == "st") "student" else "gaussian"
  for (cn in names(cases)) {
    nd <- cases[[cn]](d[d$g == "1", ])
    set.seed(77)
    mb <- brms_fill(cf, nd, 1, fam, cf$nu)[1, ]
    set.seed(77)
    mf <- frm_fill(fit, nd, cf)
    ulp <- max(abs(mb - mf) / (abs(mb) * .Machine$double.eps))
    cat(sprintf("%-7s %-12s identical %-5s max ulp %.1f\n", fn, cn,
                identical(unname(mb), unname(mf)), ulp))
  }
}
cat("\n== 2. predict() itself against brms's posterior_predict, S = 1 ==\n")
for (fn in c("ar1", "arma11", "ar2")) {
  fit <- fits[[fn]]; cf <- coefs(fit)
  nd <- cases$run4(d[d$g == "1", ])
  set.seed(78)
  mb <- brms_fill(cf, nd, 1)[1, ]
  yb <- vapply(seq_along(mb), function(i) rnorm(1, mb[i], cf$sigma), 0)
  set.seed(78)
  p <- predict(fit, newdata = nd, ndraws = 1, propagate_error = FALSE,
               summary = FALSE)
  cat(fn, "predict draw identical to brms:", identical(unname(as.vector(p)),
                                                       unname(yb)),
      " max abs diff", max(abs(as.vector(p) - yb)), "\n")
}
cat("\n== 3. grouped, distributional, S = 40000 ==\n")
S <- 40000
for (fn in c("ar1", "arma11", "ar2", "st")) {
  fit <- fits[[fn]]; cf <- coefs(fit)
  fam <- if (fn == "st") "student" else "gaussian"
  nd <- d[d$g %in% c("1", "2", "3"), ]
  nd$y[nd$g == "1" & nd$t >= 5] <- NA
  nd$y[nd$g == "2" & nd$t %in% c(1, 2, 6)] <- NA
  nd$y[nd$g == "3"] <- NA
  set.seed(79)
  mb <- brms_fill(cf, nd, S, fam, cf$nu)
  set.seed(80)
  mf <- t(vapply(seq_len(S), function(s) frm_fill(fit, nd, cf),
                 numeric(nrow(nd))))
  miss_after <- which(is.na(nd$y) | c(FALSE, head(is.na(nd$y), -1)) |
                        c(FALSE, FALSE, head(is.na(nd$y), -2)))
  zmean <- (colMeans(mf) - colMeans(mb)) /
    sqrt((apply(mf, 2, var) + apply(mb, 2, var)) / S + 1e-300)
  rsd <- apply(mf, 2, sd) / pmax(apply(mb, 2, sd), 1e-300)
  ksp <- vapply(seq_len(ncol(mb)), function(j) {
    if (sd(mb[, j]) == 0) return(NA_real_)
    suppressWarnings(ks.test(mf[, j], mb[, j])$p.value)
  }, 0)
  cat(sprintf("%s: rows with spread %d; |z| of mean max %.2f; sd ratio %.4f..%.4f; KS p min %.3g\n",
              fn, sum(!is.na(ksp)), max(abs(zmean[!is.na(ksp)])),
              min(rsd[!is.na(ksp)]), max(rsd[!is.na(ksp)]),
              min(ksp, na.rm = TRUE)))
  const <- is.na(ksp)
  cat("   rows without spread identical in both:",
      identical(unname(mb[1, const]), unname(mf[1, const])), "\n")
  # fitted(): the expected-value fill against brms's mean over fills
  fe <- fitted(fit, newdata = nd)[, "Estimate"]
  z <- (fe - colMeans(mb)) / sqrt(apply(mb, 2, var) / S + 1e-300)
  cat(sprintf("   fitted Estimate vs brms mean over fills: max |z| %.2f (rows with spread)\n",
              max(abs(z[!const]))))
  if (fn == "ar1") {
    i <- which(nd$g == "1" & nd$t >= 5)
    cat("   g1 rows 5..8 brms sd of shifted mean:",
        format(apply(mb, 2, sd)[i], digits = 5), "\n")
    cat("   theory sigma*sqrt(sum a^2k):",
        format(cf$sigma * sqrt(cumsum(cf$ar^(2 * (0:3)))) *
                 c(0, 1, 1, 1) + 0, digits = 5), "\n")
  }
}
cat("\n== 4. fitted() Est.Error on filled rows (parameter only) ==\n")
fit <- fits$ar1
nd <- d[d$g %in% c("1"), ]
nd$y[nd$t >= 5] <- NA
print(fitted(fit, newdata = nd)[5:8, ])
