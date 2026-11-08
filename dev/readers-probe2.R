# Probe of items 4, 6 and 7 on one build (READERS_ARM). Reference: Monte
# Carlo over the fit's own draw distribution at the estimates, 2e6 draws,
# seed 11 for the draws; data seeds inline.
source("dev/readers-setup.R")
ns <- asNamespace("frmtmb")
q <- function(e) tryCatch(suppressWarnings(e), error = function(err) {
  cat("ERROR:", conditionMessage(err), "\n"); NULL
})
lp_of <- function(f, nd, dp) {
  as.numeric(frm_linpred(f, newdata = nd, dpar = dp, type = "link"))
}
me_slope <- function(f, dp) {
  b <- if (dp == "mu") f$estimates$beta else {
    lp <- f$frame$linpreds[[ns$linpred_key("y", dp)]]
    f$estimates[[lp$par]][lp$idx]
  }
  nm <- if (dp == "mu") names(f$estimates$beta) else NULL
  b
}
# the draw-only variance and mean fitted() integrates, from its own route
moments <- function(f, nd) {
  ns$fitted_me_moments(f, nd, NULL, NULL, NULL, FALSE)
}
latent_data <- function(seed, ng = 8, nt = 25) {
  set.seed(seed)
  d <- expand.grid(t = 1:nt, s = factor(1:ng))
  d$x <- rnorm(nrow(d))
  d$zo <- rnorm(nrow(d))
  d$sz <- 0.5
  e <- ave(rnorm(nrow(d)), d$s, FUN = function(z) {
    stats::filter(z, 0.5, "recursive")
  })
  d$e <- 0.4 * e
  d
}

cat("\n== item 4: lognormal, latent on mu, sigma ~ me(zo, sz)\n")
d <- latent_data(9201)
zl <- d$zo + rnorm(nrow(d), 0, 0.5)
d$y <- exp(0.5 + 0.4 * d$x + d$e + rnorm(nrow(d), 0, exp(-0.7 + 0.5 * zl)))
f4 <- q(frm(bf(y ~ x + ar(time = t, gr = s, cov = TRUE), sigma ~ me(zo, sz)),
            data = d, family = lognormal()))
nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1), zo = c(0, 1.5),
                 sz = 0.5)
if (!is.null(f4)) {
  print(fitted(f4, newdata = nd))
  sl <- ns$latres_row_sd(f4, "y", NULL)
  em <- lp_of(f4, nd, "mu")
  es <- lp_of(f4, nd, "sigma")
  bs <- f4$estimates$betad
  print(bs)
  bsl <- bs[length(bs)]
  set.seed(11)
  z1 <- rnorm(2e6); z2 <- rnorm(2e6)
  mc <- sapply(1:2, function(i) {
    v <- exp(em[i] + sl * z1 + exp(2 * (es[i] + bsl * 0.5 * z2)) / 2)
    c(mean = mean(v), sd = sd(v))
  })
  cat("latent sd", sl, " me slope on sigma", bsl, "\n")
  cat("MC mean:", format(mc[1, ], digits = 6), " MC sd (draws):",
      format(mc[2, ], digits = 6), "\n")
  mm <- moments(f4, nd)
  cat("route mean:", format(mm$mean, digits = 6), " route draw sd:",
      format(sqrt(mm$var), digits = 6), "\n")
}

cat("\n== item 7b: latent + one me() on mu and sigma, lognormal\n")
d$y2 <- exp(0.5 + 0.4 * d$x + 0.3 * zl + d$e +
            rnorm(nrow(d), 0, exp(-0.7 + 0.5 * zl)))
f7b <- q(frm(bf(y2 ~ x + me(zo, sz) + ar(time = t, gr = s, cov = TRUE),
                sigma ~ me(zo, sz)), data = d, family = lognormal()))
if (!is.null(f7b)) {
  print(fitted(f7b, newdata = nd))
  sl <- ns$latres_row_sd(f7b, "y2", NULL)
  b <- f7b$estimates$beta
  bm <- b[length(b)]
  bs <- f7b$estimates$betad
  bsl <- bs[length(bs)]
  em <- as.numeric(frm_linpred(f7b, newdata = nd, dpar = "mu", type = "link"))
  es <- as.numeric(frm_linpred(f7b, newdata = nd, dpar = "sigma",
                               type = "link"))
  set.seed(11)
  z1 <- rnorm(2e6); z2 <- rnorm(2e6)
  mc <- sapply(1:2, function(i) {
    v <- exp(em[i] + bm * 0.5 * z2 + sl * z1 +
             exp(2 * (es[i] + bsl * 0.5 * z2)) / 2)
    c(mean = mean(v), sd = sd(v))
  })
  cat("latent sd", sl, " me slopes", bm, bsl, "\n")
  cat("MC mean:", format(mc[1, ], digits = 6), " MC sd (draws):",
      format(mc[2, ], digits = 6), "\n")
  mm <- moments(f7b, nd)
  cat("route mean:", format(mm$mean, digits = 6), " route draw sd:",
      format(sqrt(mm$var), digits = 6), "\n")
}

cat("\n== item 7b: latent + one me() on mu and hu, hurdle_poisson\n")
set.seed(9202)
lam <- exp(0.4 + 0.3 * d$x + 0.4 * zl + d$e)
cnt <- rpois(nrow(d), lam)
while (any(cnt == 0)) {
  i <- cnt == 0
  cnt[i] <- rpois(sum(i), lam[i])
}
d$yh <- ifelse(runif(nrow(d)) < plogis(-0.5 + 0.8 * zl), 0L, cnt)
fh <- q(frm(bf(yh ~ x + me(zo, sz) + ar(time = t, gr = s, cov = TRUE),
               hu ~ me(zo, sz)), data = d, family = hurdle_poisson()))
if (!is.null(fh)) {
  print(fitted(fh, newdata = nd))
  sl <- ns$latres_row_sd(fh, "yh", NULL)
  b <- fh$estimates$beta
  bm <- b[length(b)]
  bh <- fh$estimates$betad
  bhu <- bh[length(bh)]
  em <- as.numeric(frm_linpred(fh, newdata = nd, dpar = "mu", type = "link"))
  eh <- as.numeric(frm_linpred(fh, newdata = nd, dpar = "hu", type = "link"))
  set.seed(11)
  z1 <- rnorm(2e6); z2 <- rnorm(2e6)
  mc <- sapply(1:2, function(i) {
    l <- exp(em[i] + bm * 0.5 * z2 + sl * z1)
    v <- (1 - plogis(eh[i] + bhu * 0.5 * z2)) * l / (1 - exp(-l))
    c(mean = mean(v), sd = sd(v))
  })
  cat("latent sd", sl, " me slopes", bm, bhu, "\n")
  cat("MC mean:", format(mc[1, ], digits = 6), " MC sd (draws):",
      format(mc[2, ], digits = 6), "\n")
  mm <- moments(fh, nd)
  cat("route mean:", format(mm$mean, digits = 6), " route draw sd:",
      format(sqrt(mm$var), digits = 6), "\n")
}

cat("\n== item 7a: lognormal, me() on mu and sigma, sx .8, no latent\n")
set.seed(9203)
n <- 400
d7 <- data.frame(zo = rnorm(n), sz = 0.8)
zl <- d7$zo + rnorm(n, 0, 0.8)
d7$y <- exp(0.5 + 0.4 * zl + rnorm(n, 0, exp(-0.5 + 0.4 * zl)))
f7 <- q(frm(bf(y ~ me(zo, sz), sigma ~ me(zo, sz)), data = d7,
            family = lognormal()))
nd7 <- data.frame(zo = c(0, 1.5), sz = 0.8)
if (!is.null(f7)) {
  print(fitted(f7, newdata = nd7))
  b <- f7$estimates$beta
  bs <- f7$estimates$betad
  em <- as.numeric(frm_linpred(f7, newdata = nd7, dpar = "mu", type = "link"))
  es <- as.numeric(frm_linpred(f7, newdata = nd7, dpar = "sigma",
                               type = "link"))
  cat("me slopes", b[2], bs[2], "\n")
  set.seed(11)
  z <- rnorm(2e6)
  for (i in 1:2) {
    v <- exp(em[i] + b[2] * 0.8 * z + exp(2 * (es[i] + bs[2] * 0.8 * z)) / 2)
    cat("row", i, "MC mean", format(mean(v), digits = 6), " MC sd",
        format(sd(v), digits = 6), " max", format(max(v), digits = 4),
        " MC sd of 4 quarters:", format(sapply(split(v, rep(1:4, each = 5e5)),
                                               sd), digits = 4), "\n")
  }
  mm <- moments(f7, nd7)
  cat("route mean:", format(mm$mean, digits = 6), " route draw sd:",
      format(sqrt(mm$var), digits = 6), "\n")
}
