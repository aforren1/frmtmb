# Item 4: |ID| with by = under mm(): a path the key opens. Hand-mapped log_prob
# against frmtmb's joint at frmtmb's values. Seed 9830. Plain mm(by =) first as
# the control, then one key across mu and sigma.
source("dev/review-refuse-lib.R")
set.seed(9830)
n <- 200
lev <- paste0("L", sprintf("%02d", 1:12))
g1 <- sample(lev, n, TRUE); g2 <- sample(lev, n, TRUE)
d <- data.frame(g1 = g1, g2 = g2, x = rnorm(n))
d$f1 <- factor(ifelse(d$g1 %in% lev[1:6], "a", "b"))
d$f2 <- factor(ifelse(d$g2 %in% lev[1:6], "a", "b"))
re <- rnorm(12, 0, c(rep(0.4, 6), rep(1.2, 6))); names(re) <- lev
res <- rnorm(12, 0, 0.3); names(res) <- lev
d$y <- 1 + 0.5 * d$x + 0.5 * (re[d$g1] + re[d$g2]) + rnorm(n, 0, exp(-0.3 + 0.5 * (res[d$g1] + res[d$g2])))
txt <- c(ctl = "bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2))))",
         id = "bf(y ~ x + (1 | p | mm(g1, g2, by = cbind(f1, f2))), sigma ~ 1 + (1 | p | mm(g1, g2, by = cbind(f1, f2))))")
for (nm in names(txt)) {
  fb <- eval(parse(text = txt[[nm]]), envir = asNamespace("brms"))
  ff <- eval(parse(text = txt[[nm]]), envir = asNamespace("frmtmb"))
  pri <- E$brms_flat_prior(fb, data = d, family = stats::gaussian())
  code <- brms::make_stancode(fb, data = d, prior = pri)
  sdat <- E$brms_standata(fb, data = d, prior = pri)
  writeLines(grep("sd_1|z_1|L_1|r_1|Jby|Nby|Lcor", strsplit(code, "\n")[[1]], value = TRUE))
  pr("Jby_1 head: %s; Nby_1 %s", paste(head(sdat$Jby_1, 13), collapse = ""), sdat$Nby_1)
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
  fit <- frm(ff, data = d)
  bks <- fit$frame$re_blocks
  print(sapply(bks, function(b) c(cs = b$covstruct, dim = b$dim, nl = b$n_levels %||% NA, grp = b$group_name)))
  est <- frmtmb:::coef_b(fit); th <- fit$estimates$theta
  # one block per by level in frmtmb
  cf <- fixef(fit)
  M <- bks[[1]]$dim
  byl <- c("a", "b")
  sds <- matrix(NA_real_, M, 2)
  Ls <- array(0, c(2, M, M))
  z <- matrix(0, M, 12)
  for (k in seq_along(bks)) {
    bk <- bks[[k]]
    vc <- VarCorr(fit)[[bk$group_name]]
    b <- matrix(est[bk$c_idx], nrow = M)
    lvs <- bk$levels
    th_k <- th[bk$theta_idx]
    Sig <- bk$cov_fun %||% NULL
    s_k <- exp(th_k[seq_len(M)])
    Lk <- diag(M)
    if (M == 2L) {
      # frmtmb's us chart, by its documented map
      S <- frmtmb::us_sigma(th_k, M)
      s_k <- sqrt(diag(S))
      Lk <- t(chol(cov2cor(S)))
    }
    bylev <- k
    sds[, bylev] <- s_k
    Ls[bylev, , ] <- Lk
    z[, match(lvs, lev)] <- solve(Lk, b / s_k)
  }
  pars <- list(b = array(cf["x", 1], 1L), Intercept = cf["Intercept", 1] + cf["x", 1] * mean(d$x),
               z_1 = z)
  sn <- grep("sd_1", strsplit(code, "\n")[[1]], value = TRUE)[1]
  pars$sd_1 <- if (M == 1L && grepl("matrix", sn)) sds else sds
  if (M == 2L) pars$L_1 <- aperm(array(Ls, c(2, M, M)), c(1, 2, 3))
  if (nm == "id") pars$Intercept_sigma <- cf["sigma_Intercept", 1]
  if (nm == "ctl") pars$sigma <- exp(fit$estimates$betad[[1]])
  up <- tr(rstan::unconstrain_pars(sf, pars))
  if (is_err(up)) { pr("%s unconstrain: %s", nm, msg(up, 200)); print(sf@par_dims); next }
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
  # frmtmb's joint: the MVN of b, brms's std_normal(z): log det per level
  j <- -fit$obj$env$f(fit$obj$env$last.par.best)
  chg <- 0
  for (k in seq_along(bks)) {
    S <- diag(sds[, k], M) %*% matrix(Ls[k, , ], M, M) %*% t(diag(sds[, k], M) %*% matrix(Ls[k, , ], M, M))
    chg <- chg + bks[[k]]$n_levels * 0.5 * as.numeric(determinant(S)$modulus)
  }
  g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  p3 <- pars; p3$z_1 <- p3$z_1 * 0.9
  gz <- g[rstan::unconstrain_pars(sf, p3) != up]
  pr("%-4s brms lp %.6f, frmtmb joint %.6f, change of variables %.6f; difference %.3e; inner grad max %.3e",
     nm, lp, j, chg, lp - j - chg, max(abs(gz)))
}
