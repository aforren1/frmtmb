# The route decision, measured: frm()'s ML sigma under a cov = TRUE
# ar(), ma(), arma() (order 1), ar(p = 2), and ar() with sigma ~ x.
#
# Arms:
#   before  rellib-r14 (frmtmb 0.76.0, sigma the marginal sd). Prints the
#           fit and "route B": the marginal estimate carried to the
#           innovation sd, sigma * c(theta), with a delta-method SE from
#           vcov(full = TRUE).
#   after   the lane library (sigma the innovation sd, "route A").
#   ref     an independent ML fit on the innovation scale: a dense
#           multivariate normal written from brms's covariance
#           (cholesky_cor_*() of brms 2.23.0, read in
#           dev/sigmainn-out/brms-code.txt), optim() then a numDeriv
#           Hessian. No frmtmb code.
#   summary reads the three RDS files and prints the comparison.
# Data (seed 20261027): ar, the AR(1) 0.9 data of dev/acprior-agree.R;
# ma and arma as in dev/acprior-innov.R; ar2, an AR(2) with phi 0.5,
# 0.3 in 15 groups of 8; sigx, the ar data with log sigma_innov =
# -0.8 + 0.3 x2 (x2 a second covariate).
# Usage: Rscript dev/sigmainn-route.R before|after|ref|summary
# Log: dev/sigmainn-out/route-<arm>.txt, RDS route-<arm>.rds
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
out <- file.path(wt, "dev", "sigmainn-out")
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14", ulib))
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
sim_arma <- function(ng, k, ar, ma, seed) {
  set.seed(seed)
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                  x = rnorm(ng * k))
  e <- unlist(lapply(seq_len(ng), function(i) {
    as.numeric(stats::arima.sim(list(ar = if (any(ar != 0)) ar,
                                      ma = if (ma != 0) ma), n = k))
  }))
  d$y <- 1 + 0.5 * d$x + e
  d
}
mk_sigx <- function() {
  set.seed(20261027)
  ng <- 20; k <- 6
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                  x = rnorm(ng * k), x2 = rnorm(ng * k))
  si <- exp(-0.8 + 0.3 * d$x2)
  e <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    # an AR(1) 0.7 innovation process with a row-varying innovation sd,
    # brms's model: e = D z, z with correlation R / (1 - ar^2)
    e[i] <- if (d$t[i] == 1) rnorm(1, 0, 1 / sqrt(1 - 0.49)) else
      0.7 * e[i - 1] + rnorm(1)
  }
  d$y <- 1 + 0.5 * d$x + si * e
  d
}
dat <- list(ar = make_data()$ar,
            ma = sim_arma(10, 6, 0, 0.6, 20261027),
            arma = sim_arma(12, 6, 0.7, 0.3, 20261027),
            ar2 = sim_arma(15, 8, c(0.5, 0.3), 0, 20261027),
            sigx = mk_sigx())
fo <- list(ar = y ~ x + ar(time = t, gr = g, cov = TRUE),
           ma = y ~ x + ma(time = t, gr = g, cov = TRUE),
           arma = y ~ x + arma(time = t, gr = g, cov = TRUE),
           ar2 = y ~ x + ar(time = t, gr = g, p = 2, cov = TRUE),
           sigx = NULL)
pq <- list(ar = c(1, 0), ma = c(0, 1), arma = c(1, 1), ar2 = c(2, 0),
           sigx = c(1, 0))

# brms's covariance over sigma_innov^2 (its cholesky_cor_*() before the
# Cholesky), and for AR(2) the autocovariance of unit innovations from
# the psi weights, summed to 5000 terms
brms_cov <- function(ar, ma, n) {
  if (length(ar) == 2L) {
    psi <- c(1, stats::ARMAtoMA(ar = ar, lag.max = 5000))
    acv <- vapply(0:(n - 1), function(h) {
      sum(psi[1:(length(psi) - h)] * psi[(1 + h):length(psi)])
    }, 0)
    return(stats::toeplitz(acv))
  }
  ar <- if (length(ar)) ar else 0
  ma <- if (length(ma)) ma else 0
  if (ma == 0) return(ar^abs(outer(1:n, 1:n, "-")) / (1 - ar^2))
  if (ar == 0) {
    m <- diag(1 + ma^2, n)
    m[abs(row(m) - col(m)) == 1] <- ma
    return(m)
  }
  g1 <- (1 + ar * ma) * (ar + ma)
  m <- diag(1 + ma^2 + 2 * ar * ma, n)
  for (i in 2:n) for (j in 1:(i - 1)) {
    m[i, j] <- m[j, i] <- g1 * ar^(i - j - 1)
  }
  m / (1 - ar^2)
}
# the marginal sd over the innovation sd, from the same covariance
cfac <- function(ar, ma) 1 / sqrt(brms_cov(ar, ma, 2)[1, 1])

if (arm %in% c("before", "after")) {
  suppressPackageStartupMessages(library(frmtmb))
  cat("frmtmb", format(packageVersion("frmtmb")), "from",
      dirname(find.package("frmtmb")), "\n")
  res <- list()
  for (m in names(dat)) {
    d <- dat[[m]]
    fit <- if (m == "sigx") {
      frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE), sigma ~ x2), data = d)
    } else frm(fo[[m]], data = d)
    V <- vcov(fit, full = TRUE)
    est <- fit$obj$env$last.par.best
    nat <- frmtmb:::autocor_natural(fit$estimates$thetaac,
                                    fit$frame$autocor[[1]])
    ar <- nat[grepl("^ar", names(nat))]
    ma <- nat[grepl("^ma", names(nat))]
    cc <- cfac(ar, ma)
    th <- fit$estimates$thetaac
    lsrow <- grep("^sigma_\\(Intercept\\)$", rownames(V))
    throws <- grep("^thetaac", rownames(V))
    # d log c / d theta by central differences of the reference c
    dlc <- vapply(seq_along(th), function(j) {
      h <- 1e-6 * max(1, abs(th[j]))
      tp <- th; tp[j] <- tp[j] + h
      tm <- th; tm[j] <- tm[j] - h
      np <- frmtmb:::autocor_natural(tp, fit$frame$autocor[[1]])
      nm <- frmtmb:::autocor_natural(tm, fit$frame$autocor[[1]])
      (log(cfac(np[grepl("^ar", names(np))], np[grepl("^ma", names(np))])) -
         log(cfac(nm[grepl("^ar", names(nm))], nm[grepl("^ma", names(nm))]))) /
        (2 * h)
    }, 0)
    J <- numeric(nrow(V)); J[lsrow] <- 1; J[throws] <- dlc
    ls <- fit$estimates$betad[["sigma_(Intercept)"]]
    r <- list(
      logLik = as.numeric(logLik(fit)),
      beta = fit$estimates$beta, betad = fit$estimates$betad,
      thetaac = th, nat = nat, c = cc,
      sigma = if (m != "sigx") sigma(fit) else NA,
      se_logsigma = sqrt(V[lsrow, lsrow]),
      se_beta = sqrt(diag(V))[seq_along(fit$estimates$beta)],
      se_theta = sqrt(diag(V))[throws],
      cor_ls_theta = V[lsrow, throws[1]] /
        sqrt(V[lsrow, lsrow] * V[throws[1], throws[1]]),
      # route B: the marginal estimate carried to the innovation sd
      routeB_logsigma = ls + log(cc),
      routeB_se = sqrt(drop(t(J) %*% V %*% J)),
      routeB_cor = drop(t(J) %*% V[, throws[1]]) /
        sqrt(drop(t(J) %*% V %*% J) * V[throws[1], throws[1]]),
      iter = fit$opt$iterations %||% NA,
      evals = paste(fit$opt$evaluations, collapse = "/"),
      summary_sigma = if (m != "sigx") summary(fit)$spec_pars else NULL,
      varcorr = if (m != "sigx") VarCorr(fit)$residual__$sd else NULL,
      slope = if (m == "sigx") fit$estimates$betad[["sigma_x2"]] else NA)
    res[[m]] <- r
    cat(sprintf(paste0("%-5s logLik %.15g | sigma_(Intercept) %.12f se %.6f",
                       " | c %.10f | routeB %.12f se %.6f | cor(ls,th) ",
                       "%.3f routeB %.3f | iter %s evals %s\n"),
                m, r$logLik, ls, r$se_logsigma, cc, r$routeB_logsigma,
                r$routeB_se, r$cor_ls_theta, r$routeB_cor, r$iter,
                r$evals))
  }
  saveRDS(res, file.path(out, paste0("route-", arm, ".rds")))
} else if (identical(arm, "ref")) {
  suppressPackageStartupMessages(library(numDeriv))
  res <- list()
  for (m in names(dat)) {
    d <- dat[[m]]
    p <- pq[[m]][1]; q <- pq[[m]][2]
    X <- cbind(1, d$x)
    Xs <- if (m == "sigx") cbind(1, d$x2) else matrix(1, nrow(d), 1)
    gs <- split(seq_len(nrow(d)), d$g)
    unpack <- function(par) {
      b <- par[1:2]
      ls <- par[2 + seq_len(ncol(Xs))]
      k <- 2 + ncol(Xs)
      # stationary AR(2) through its partial autocorrelations, the
      # Durbin-Levinson map; ar(1) and ma(1) through tanh
      if (p == 2) {
        r <- tanh(par[k + 1:2])
        ar <- c(r[1] * (1 - r[2]), r[2])
        ma <- numeric(0)
      } else {
        ar <- if (p) tanh(par[k + 1]) else numeric(0)
        ma <- if (q) tanh(par[k + p + 1]) else numeric(0)
      }
      list(b = b, ls = ls, ar = ar, ma = ma)
    }
    nll <- function(par) {
      u <- unpack(par)
      mu <- drop(X %*% u$b)
      si <- exp(drop(Xs %*% u$ls))
      out <- 0
      for (rows in gs) {
        k <- length(rows)
        S <- brms_cov(u$ar, u$ma, k) * outer(si[rows], si[rows])
        L <- tryCatch(chol(S), error = function(e) NULL)
        # a step of optim() past the edge, where tanh() rounds to 1
        if (is.null(L) || anyNA(L)) return(1e10)
        z <- backsolve(L, y <- d$y[rows] - mu[rows], transpose = TRUE)
        out <- out + sum(log(diag(L))) + 0.5 * sum(z^2) +
          0.5 * k * log(2 * pi)
      }
      out
    }
    st <- c(stats::coef(lm(y ~ x, d)), log(stats::sd(d$y)),
            if (ncol(Xs) > 1) 0, rep(0.2, p + q))
    o <- stats::optim(st, nll, method = "BFGS",
                      control = list(maxit = 5000, reltol = 1e-15))
    # Newton steps on the numDeriv Hessian, to the optimizer's precision
    par <- o$par
    for (it in 1:20) {
      g <- numDeriv::grad(nll, par)
      H <- numDeriv::hessian(nll, par)
      step <- solve(H, g)
      par <- par - step
      if (max(abs(step)) < 1e-12) break
    }
    H <- numDeriv::hessian(nll, par)
    V <- solve(H)
    u <- unpack(par)
    r <- list(logLik = -nll(par), b = u$b, logsigma = u$ls,
              se_logsigma = sqrt(diag(V))[2 + seq_along(u$ls)],
              ar = u$ar, ma = u$ma,
              se_b = sqrt(diag(V))[1:2], grad = max(abs(numDeriv::grad(nll, par))))
    res[[m]] <- r
    cat(sprintf("%-5s logLik %.15g | log sigma_innov %s se %s | ar %s ma %s | max|grad| %.2e\n",
                m, r$logLik, paste(sprintf("%.12f", u$ls), collapse = ","),
                paste(sprintf("%.6f", r$se_logsigma), collapse = ","),
                paste(sprintf("%.8f", u$ar), collapse = ","),
                paste(sprintf("%.8f", u$ma), collapse = ","), r$grad))
  }
  saveRDS(res, file.path(out, "route-ref.rds"))
} else {
  B <- readRDS(file.path(out, "route-before.rds"))
  A <- if (file.exists(file.path(out, "route-after.rds")))
    readRDS(file.path(out, "route-after.rds"))
  R <- readRDS(file.path(out, "route-ref.rds"))
  rel <- function(a, b) abs(a - b) / max(abs(b), 1e-300)
  for (m in names(B)) {
    cat("\n=====", m, "\n")
    b <- B[[m]]; r <- R[[m]]
    cat(sprintf("logLik: before %.15g | ref %.15g | rel %.2e\n",
                b$logLik, r$logLik, rel(b$logLik, r$logLik)))
    if (!is.null(A)) {
      a <- A[[m]]
      cat(sprintf("logLik: after  %.15g | rel to before %.2e\n",
                  a$logLik, rel(a$logLik, b$logLik)))
    }
    cat(sprintf("log sigma_innov (intercept): route B %.12f | ref %.12f | rel %.2e\n",
                b$routeB_logsigma, r$logsigma[1],
                rel(b$routeB_logsigma, r$logsigma[1])))
    if (!is.null(A)) {
      cat(sprintf("log sigma_innov (intercept): route A %.12f | rel to ref %.2e | rel to B %.2e\n",
                  A[[m]]$betad[["sigma_(Intercept)"]], rel(A[[m]]$betad[["sigma_(Intercept)"]], r$logsigma[1]),
                  rel(A[[m]]$betad[["sigma_(Intercept)"]], b$routeB_logsigma)))
    }
    cat(sprintf("SE log sigma_innov: route B %.6f | ref %.6f | ratio %.6f\n",
                b$routeB_se, r$se_logsigma[1], b$routeB_se / r$se_logsigma[1]))
    if (!is.null(A)) {
      cat(sprintf("SE log sigma_innov: route A %.6f | ratio to ref %.6f | to B %.6f\n",
                  A[[m]]$se_logsigma, A[[m]]$se_logsigma / r$se_logsigma[1],
                  A[[m]]$se_logsigma / b$routeB_se))
      cat(sprintf("cor(log sigma, theta_1): before (marginal) %.3f | route B %.3f | after %.3f\n",
                  b$cor_ls_theta, b$routeB_cor, A[[m]]$cor_ls_theta))
      cat(sprintf("optimizer iterations/evals: before %s %s | after %s %s\n",
                  b$iter, b$evals, A[[m]]$iter, A[[m]]$evals))
      cat(sprintf("beta: max rel after vs before %.2e | after vs ref %.2e\n",
                  max(abs(A[[m]]$beta - b$beta) / abs(b$beta)),
                  max(abs(A[[m]]$beta - r$b) / abs(r$b))))
      if (m == "sigx") {
        cat(sprintf("sigma slope: before %.10f after %.10f ref %.10f\n",
                    b$slope, A[[m]]$slope, r$logsigma[2]))
      } else {
        cat(sprintf("sigma(fit): before %.10f after %.10f ref %.10f; summary %.4f -> %.4f; VarCorr %.4f -> %.4f\n",
                    b$sigma, A[[m]]$sigma, exp(r$logsigma[1]),
                    b$summary_sigma$Estimate, A[[m]]$summary_sigma$Estimate,
                    b$varcorr[1], A[[m]]$varcorr[1]))
      }
    }
    cat(sprintf("natural ar/ma: before %s | ref %s\n",
                paste(sprintf("%.8f", b$nat), collapse = ","),
                paste(sprintf("%.8f", c(r$ar, r$ma)), collapse = ",")))
  }
}
