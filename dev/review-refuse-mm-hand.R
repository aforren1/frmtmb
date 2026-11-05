# Item 4: the translator stops on a weighted mm(), so the lane has no L1 on
# weights under a key or under cov. By hand: brms's log_prob at frmtmb's values.
# Seed 9840.
source("dev/review-refuse-lib.R")
set.seed(9840)
n <- 150
lev <- paste0("L", sprintf("%02d", 1:12))
g1 <- sample(lev[1:10], n, TRUE, prob = c(0.35, rep(0.65 / 9, 9)))
g2 <- sample(lev[2:12], n, TRUE, prob = c(rep(0.6 / 9, 9), 0.2, 0.2))
d <- data.frame(g1 = g1, g2 = g2, x = rnorm(n), w1 = runif(n, 0.2, 1), w2 = runif(n, 0.2, 1),
                x1 = rnorm(n), x2 = rnorm(n))
re <- rnorm(12, 0, 0.8); names(re) <- lev; res <- rnorm(12, 0, 0.3); names(res) <- lev
d$y <- 1 + 0.5 * d$x + (d$w1 * re[d$g1] + d$w2 * re[d$g2]) / (d$w1 + d$w2) +
  rnorm(n, 0, exp(-0.2 + (d$w1 * res[d$g1] + d$w2 * res[d$g2]) / (d$w1 + d$w2)))
A <- exp(-abs(outer(1:12, 1:12, "-")) / 3); dimnames(A) <- list(lev, lev)
L <- t(chol(A))
cases <- list(
  cov_w = list(txt = "bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2), cov = A)))", M = 1L, cov = TRUE, sig = FALSE),
  cov_w2 = list(txt = "bf(y ~ x + (1 + x1 | mm(g1, g2, weights = cbind(w1, w2), cov = A)))", M = 2L, cov = TRUE, sig = FALSE),
  id_w = list(txt = "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))", M = 2L, cov = FALSE, sig = TRUE),
  id_cov_w = list(txt = "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2), cov = A)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2), cov = A)))", M = 2L, cov = TRUE, sig = TRUE))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fb <- eval(parse(text = cs$txt), envir = asNamespace("brms"))
  ff <- eval(parse(text = cs$txt), envir = asNamespace("frmtmb"))
  pri <- E$brms_flat_prior(fb, data = d, family = stats::gaussian(), data2 = list(A = A))
  code <- brms::make_stancode(fb, data = d, prior = pri, data2 = list(A = A))
  sdat <- E$brms_standata(fb, data = d, prior = pri, data2 = list(A = A))
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
  fit <- tr(frm(ff, data = d, data2 = list(A = A)))
  if (is_err(fit)) { pr("%s %s", nm, msg(fit, 150)); next }
  bk <- fit$frame$re_blocks[[1]]
  est <- frmtmb:::coef_b(fit)
  cf <- fixef(fit)
  M <- cs$M
  B <- matrix(est[bk$c_idx], nrow = M)
  lv <- bk$levels
  th <- fit$estimates$theta[bk$theta_idx]
  if (M == 1L) {
    s <- exp(th[1]); Lc <- matrix(1, 1, 1)
  } else if (identical(bk$covstruct, "gr_cov") || TRUE) {
    S <- frmtmb::us_sigma(th, M)
    s <- sqrt(diag(S)); Lc <- t(chol(cov2cor(S)))
  }
  # brms: r = (D Lc z (L_A)')', one row per level
  LA <- if (cs$cov) t(chol(A[lv, lv])) else diag(length(lv))
  z <- solve(Lc, B / s) %*% t(solve(LA))
  pars <- list(b = array(cf["x", 1], 1L), Intercept = cf["Intercept", 1] + cf["x", 1] * mean(d$x),
               sd_1 = array(s, M))
  pars$z_1 <- if (grepl("array\\[M_1\\] vector\\[N_1\\] z_1", code)) lapply(seq_len(M), function(i) z[i, ]) else z
  if (M == 2L) pars$L_1 <- Lc
  if (cs$sig) pars$Intercept_sigma <- cf["sigma_Intercept", 1] else pars$sigma <- exp(fit$estimates$betad[[1]])
  up <- tr(rstan::unconstrain_pars(sf, pars))
  if (is_err(up)) { pr("%s: %s", nm, msg(up, 220)); print(names(sf@par_dims)); next }
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
  j <- -fit$obj$env$f(fit$obj$env$last.par.best)
  Sg <- diag(s, M) %*% Lc %*% t(Lc) %*% diag(s, M)
  chg <- 0.5 * (length(lv) * as.numeric(determinant(Sg)$modulus) + M * if (cs$cov) as.numeric(determinant(A[lv, lv])$modulus) else 0)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  p3 <- pars
  p3$z_1 <- if (is.list(pars$z_1)) lapply(pars$z_1, function(v) v * 0.9) else pars$z_1 * 0.9
  gz <- g[rstan::unconstrain_pars(sf, p3) != up]
  pr("%-9s brms lp %.6f; frmtmb joint %.6f; change of variables %.6f; difference %.3e; inner grad max %.3e", nm, lp, j, chg, lp - j - chg, max(abs(gz)))
}
