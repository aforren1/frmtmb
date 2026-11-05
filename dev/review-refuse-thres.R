# Item 6: thres() x mvbf with a key across the responses, by hand, because the
# translator stops on an ordinal response in a multivariate model. Seed 10003.
source("dev/review-refuse-lib.R")
set.seed(10003); n <- 180; x <- rnorm(n); g <- factor(sample(1:10, n, TRUE))
u <- rnorm(10, 0, 0.6)
o <- cut(0.7 * x + u[g] + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
# a correlation away from its edge, or the Cholesky factor reads rounding
v <- 0.4 * u + rnorm(10, 0, 0.6)
d <- data.frame(x = x, o = o, w = rnorm(n, 1 + 0.3 * x + v[g]), g = g)
fb <- brms::bf(o | thres(3) ~ x + (1 | p | g), family = brms::cumulative()) + brms::bf(w ~ x + (1 | p | g)) + brms::set_rescor(FALSE)
pri <- E$brms_flat_prior(fb, data = d, family = NULL)
code <- brms::make_stancode(fb, data = d, prior = pri)
sdat <- E$brms_standata(fb, data = d, prior = pri)
sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
print(names(sf@par_dims)); print(unlist(lapply(sf@par_dims, function(v) paste(v, collapse = "x"))))
fit <- frm(mvbf(bf(o | thres(3) ~ x + (1 | p | g), family = cumulative()), bf(w ~ x + (1 | p | g)), rescor = FALSE), data = d)
print(fixef(fit)[, 1])
print(sapply(fit$frame$re_blocks, function(b) c(b$covstruct, b$dim, b$group_name, paste(b$cnms, collapse = ","))))
cf <- fixef(fit)[, 1]
est <- frmtmb:::coef_b(fit)
bk <- fit$frame$re_blocks[[1]]
B <- matrix(est[bk$c_idx], nrow = bk$dim)
S <- frmtmb::us_sigma(fit$estimates$theta[bk$theta_idx], 2L)
s <- sqrt(diag(S)); Lc <- t(chol(cov2cor(S)))
z <- solve(Lc, B / s)
nm <- names(cf)
pr("names: %s", paste(nm, collapse = " "))
ths <- cf[grep("Intercept.*\\[", nm)]
if (!length(ths)) ths <- cf[grep("^o_Intercept", nm)]
pr("thresholds: %s", paste(round(ths, 3), collapse = " "))
xc <- mean(d$x)
bo <- cf[["o_x"]]; bw <- cf[["w_x"]]
pars <- list(b_o = array(bo, 1L), Intercept_o = array(ths - bo * xc, length(ths)),
             b_w = array(bw, 1L), Intercept_w = cf[["w_Intercept"]] + bw * xc,
             sigma_w = exp(fit$estimates$betad[[length(fit$estimates$betad)]]),
             sd_1 = array(s, 2L), L_1 = Lc, z_1 = z)
up <- tr(rstan::unconstrain_pars(sf, pars))
if (is_err(up)) stop(msg(up, 300))
lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
j <- -fit$obj$env$f(fit$obj$env$last.par.best)
chg <- 0.5 * nlevels(g) * as.numeric(determinant(S)$modulus)
g2 <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
p3 <- pars; p3$z_1 <- p3$z_1 * 0.9
gz <- g2[rstan::unconstrain_pars(sf, p3) != up]
pr("thres + key: brms lp %.6f; frmtmb joint %.6f; change of variables %.6f; difference %.3e; inner grad max %.3e; ordered threshold Jacobian not in play under adjust_transform = FALSE",
   lp, j, chg, lp - j - chg, max(abs(gz)))
