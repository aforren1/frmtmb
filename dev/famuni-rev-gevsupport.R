# Reviewer probe: does a GEV fit end inside the support, is logLik the
# true log-likelihood, and does it agree with evd::fgev?
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
true_ll <- function(y, mu, sg, xi) {
  sum(mapply(function(y, m, s, k) evd::dgev(y, m, s, k, log = TRUE),
             y, mu, sg, xi))
}
tmin <- function(fit, y) {
  mu <- frm_linpred(fit, dpar = "mu", type = "response")
  sg <- frm_linpred(fit, dpar = "sigma", type = "response")
  xi <- frm_linpred(fit, dpar = "xi", type = "response")
  t <- 1 + xi * (y - mu) / sg
  list(t = t, mu = mu, sg = sg, xi = xi)
}
wcatch <- function(expr) {
  ws <- character(0); ms <- character(0)
  v <- withCallingHandlers(expr,
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    },
    message = function(m) {
      ms <<- c(ms, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  list(v = v, w = ws, m = ms)
}

cat("\n=== A. no covariates, against evd::fgev\n")
for (cfg in list(c(-0.45, 1), c(-0.7, 2), c(-0.9, 3), c(-0.3, 4),
                 c(0.3, 5), c(0.8, 6))) {
  set.seed(400 + cfg[2])
  n <- 150
  y <- rg(n, 2, 1.5, cfg[1])
  r <- wcatch(frm(y ~ 1, family = gen_extreme_value(),
                  data = data.frame(y = y)))
  fit <- r$v
  tm <- tmin(fit, y)
  ev <- suppressWarnings(evd::fgev(y))
  ll_true <- true_ll(y, tm$mu, tm$sg, tm$xi)
  cat(sprintf(
    "xi_true %5.2f | frm: conv %d logLik %.6f xi %.5f min t %.3e | true ll at frm %.6f | fgev: ll %.6f xi %.5f (se %s) conv %d\n",
    cfg[1], fit$opt$convergence, as.numeric(logLik(fit)), tm$xi[1],
    min(tm$t), ll_true, -ev$deviance / 2, ev$estimate[["shape"]],
    format(signif(ev$std.err[["shape"]], 4)), ev$convergence == "successful"))
  sm <- summary(fit)
  fx <- fixef(fit)
  cat("   fixef rows:", paste(rownames(fx), collapse = " "), "\n")
  cat("   fixef xi row (link scale):",
      paste(signif(fx[grep("xi", rownames(fx)), 1:2], 5), collapse = " "),
      " evd se mu/sigma/xi:",
      paste(signif(ev$std.err, 4), collapse = " "), "\n")
  if (length(r$w)) cat("   WARN:", substr(r$w, 1, 200), sep = "\n   ")
  if (length(r$m)) cat("   MSG:", substr(r$m, 1, 200), sep = "\n   ")
  cat("\n")
}

cat("\n=== B. covariates and a random intercept, xi near -0.5 and -0.8\n")
for (cfg in list(c(-0.5, 11), c(-0.8, 12), c(-0.95, 13), c(0.5, 14))) {
  set.seed(500 + cfg[2])
  n <- 400
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = gl(20, n / 20))
  b <- stats::rnorm(20, 0, 0.4)[d$g]
  d$y <- rg(n, 1 + 0.5 * d$x + b, exp(-0.2 + 0.3 * d$z), cfg[1])
  r <- wcatch(frm(bf(y ~ x + (1 | g), sigma ~ z),
                  family = gen_extreme_value(), data = d))
  fit <- r$v
  tm <- tmin(fit, d$y)
  cat(sprintf("xi_true %5.2f | conv %d logLik %.6f xi %.5f | min t %.3e | rows t<=0: %d | rows t<1e-6: %d\n",
              cfg[1], fit$opt$convergence, as.numeric(logLik(fit)),
              tm$xi[1], min(tm$t), sum(tm$t <= 0), sum(tm$t < 1e-6)))
  print(signif(fixef(fit)[, 1:2], 5))
  if (length(r$w)) cat("   WARN:", substr(r$w, 1, 300), sep = "\n   ")
  if (length(r$m)) cat("   MSG:", substr(r$m, 1, 300), sep = "\n   ")
  # the same model without the random effect, so the conditional
  # log-likelihood is the reported one, against evd's density
  r2 <- wcatch(frm(bf(y ~ x, sigma ~ z), family = gen_extreme_value(),
                   data = d))
  tm2 <- tmin(r2$v, d$y)
  cat(sprintf("   no-RE fit: logLik %.8f, sum evd::dgev at its estimates %.8f, min t %.3e\n",
              as.numeric(logLik(r2$v)),
              true_ll(d$y, tm2$mu, tm2$sg, tm2$xi), min(tm2$t)))
  cat("\n")
}

cat("\n=== C. link_xi = 'identity', where xi < -1 is reachable\n")
for (cfg in list(c(-0.8, 21), c(-0.95, 22), c(-0.6, 23))) {
  set.seed(600 + cfg[2])
  n <- 60
  y <- rg(n, 0, 1, cfg[1])
  r <- wcatch(frm(y ~ 1, family = gen_extreme_value(link_xi = "identity"),
                  data = data.frame(y = y)))
  fit <- r$v
  tm <- tmin(fit, y)
  ev <- suppressWarnings(evd::fgev(y, std.err = FALSE))
  cat(sprintf(
    "xi_true %5.2f | identity: conv %d logLik %.4f xi %.4f min t %.3e true ll %.4f | fgev ll %.4f xi %.4f\n",
    cfg[1], fit$opt$convergence, as.numeric(logLik(fit)), tm$xi[1],
    min(tm$t), true_ll(y, tm$mu, tm$sg, tm$xi), -ev$deviance / 2,
    ev$estimate[["shape"]]))
  if (length(r$w)) cat("   WARN:", substr(r$w, 1, 300), sep = "\n   ")
  # a start beyond -1: the log-likelihood at xi = -1.5 with the upper end
  # placed just above the largest observation, from the family's density
  fam <- gen_extreme_value(link_xi = "identity")
  for (xi0 in c(-1.5, -3)) {
    for (gap in c(1e-6, 1e-12, 0, -1e-9)) {
      mu0 <- mean(y); s0 <- (max(y) + gap - mu0) * (-xi0)
      lp <- fam$lpdf(y, list(mu = mu0, sigma = s0, xi = xi0,
                             .eta_sigma = log(s0)), list())
      cat(sprintf("   xi %.1f gap %g: family logLik %.4f (max row %.3f)\n",
                  xi0, gap, sum(lp), max(lp)))
    }
  }
}

cat("\n=== D. how a scalar xi is reported\n")
set.seed(701)
y <- rg(200, 2, 1.5, -0.3)
fit <- frm(y ~ 1, family = gen_extreme_value(), data = data.frame(y = y))
print(summary(fit))
cat("fixef:\n"); print(fixef(fit))
cat("coef names:", names(coef(fit)), "\n")
cat("frm_linpred xi response [1]:",
    frm_linpred(fit, dpar = "xi", type = "response")[1], "\n")
