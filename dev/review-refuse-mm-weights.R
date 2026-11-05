# Item 4: terms with one key and DIFFERENT weights. The lane says each term
# keeps its own weights, as a claim of parity. brms's standata holds W_1_1 and
# W_1_2 per group id. Hand log_prob at chosen values under each reading, and
# with the terms swapped; frmtmb's reading at its own fit. Seed 9810.
source("dev/review-refuse-lib.R")
set.seed(9810)
n <- 120
lev <- paste0("L", sprintf("%02d", 1:8))
d <- data.frame(g1 = sample(lev[1:6], n, TRUE), g2 = sample(lev[3:8], n, TRUE),
                x = rnorm(n), w1 = runif(n, 0.2, 1), w2 = runif(n, 0.2, 1),
                v1 = runif(n, 0.2, 1), v2 = runif(n, 0.2, 1))
re <- rnorm(8, 0, 0.8); names(re) <- lev
res <- rnorm(8, 0, 0.5); names(res) <- lev
d$y <- 1 + 0.5 * d$x + (d$w1 * re[d$g1] + d$w2 * re[d$g2]) / (d$w1 + d$w2) +
  rnorm(n, 0, exp((d$v1 * res[d$g1] + d$v2 * res[d$g2]) / (d$v1 + d$v2)))
txt <- list(
  wv = "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2))))",
  vw = "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(v1, v2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
set.seed(9811)
sds <- c(0.5, 0.3); rho <- 0.4
Lc <- t(chol(matrix(c(1, rho, rho, 1), 2)))
z <- matrix(rnorm(16), 2, 8)
B <- diag(sds) %*% Lc %*% z; colnames(B) <- lev
b0 <- 0.9; b1 <- 0.4; s0 <- -0.1
wmean <- function(a1, a2, r) (a1 * r[d$g1] + a2 * r[d$g2]) / (a1 + a2)
for (nm in names(txt)) {
  fb <- eval(parse(text = txt[[nm]]), envir = asNamespace("brms"))
  pri <- E$brms_flat_prior(fb, data = d, family = stats::gaussian())
  code <- brms::make_stancode(fb, data = d, prior = pri)
  sdat <- E$brms_standata(fb, data = d, prior = pri)
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
  pars <- list(b = array(b1, 1L), Intercept = b0 + b1 * mean(d$x), Intercept_sigma = s0,
               sd_1 = array(sds, 2L), z_1 = z, L_1 = Lc)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
  zd <- sum(dnorm(z, log = TRUE))
  wa <- if (nm == "wv") list(d$w1, d$w2) else list(d$v1, d$v2)
  va <- if (nm == "wv") list(d$v1, d$v2) else list(d$w1, d$w2)
  muh <- b0 + b1 * d$x + wmean(wa[[1]], wa[[2]], B[1, ])
  own <- sum(dnorm(d$y, muh, exp(s0 + wmean(va[[1]], va[[2]], B[2, ])), log = TRUE)) + zd
  first <- sum(dnorm(d$y, muh, exp(s0 + wmean(wa[[1]], wa[[2]], B[2, ])), log = TRUE)) + zd
  pr("%s: brms log_prob %.6f; hand, each term its own weights %.6f (diff %.2e); hand, the first term's weights for both %.6f (diff %.2e)",
     nm, lp, own, lp - own, first, lp - first)
}
# frmtmb's reading
ff <- eval(parse(text = txt$wv), envir = asNamespace("frmtmb"))
fit <- frm(ff, data = d)
cf <- fixef(fit)
est <- frmtmb:::coef_b(fit)
bk <- fit$frame$re_blocks[[1]]
Bf <- matrix(est[bk$c_idx], nrow = bk$dim)
colnames(Bf) <- bk$levels %||% dimnames(ranef(fit)[[1]])[[1]]
pr("frmtmb block %s, %s rows", bk$group_name, nrow(Bf))
sg <- as.numeric(fitted(fit, dpar = "sigma")[, 1])
own <- exp(cf["sigma_Intercept", 1] + wmean(d$v1, d$v2, Bf[2, ]))
fst <- exp(cf["sigma_Intercept", 1] + wmean(d$w1, d$w2, Bf[2, ]))
pr("frmtmb fitted sigma vs hand: own weights %.3e, first term's weights %.3e", max(abs(sg - own)), max(abs(sg - fst)))
# the two models' fits differ by how much on these data
fit1 <- frm(bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2)))), data = d)
pr("logLik frm own weights %.4f; frm with brms's model (the first term's for both) %.4f", as.numeric(logLik(fit)), as.numeric(logLik(fit1)))
# and brms with weights on the first term and no weights on the second
txt2 <- "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2)))"
sd2 <- brms::make_standata(eval(parse(text = txt2), envir = asNamespace("brms")), data = d)
pr("first weighted, second not: W_1_1 head %s (a plain mm() is 0.5)", paste(round(head(sd2$W_1_1, 3), 3), collapse = " "))
txt3 <- "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))"
sd3 <- brms::make_standata(eval(parse(text = txt3), envir = asNamespace("brms")), data = d)
pr("second weighted, first not: W_1_1 head %s", paste(round(head(sd3$W_1_1, 3), 3), collapse = " "))
# mmc columns under one key: Z, with weights
txt4 <- "bf(y ~ x + (1 + mmc(x, w1) | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2)))"
sd4 <- brms::make_standata(eval(parse(text = txt4), envir = asNamespace("brms")), data = d)
print(grep("^(W|Z)_", names(sd4), value = TRUE))
