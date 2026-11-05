# Item 4: the me() draw's variance on the mean, per family route, against brms's
# draws at the same values. Seed 7971. A row is a route through fitted(): the mean
# is a link inverse, a mean function of several dpars, or category probabilities.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7971)
n <- 400
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cnt <- rpois(n, exp(eta))
d$nb <- rnbinom(n, mu = exp(eta), size = 3)
d$pos <- exp(eta + rnorm(n, 0, 0.5))
d$gam <- rgamma(n, 4, 4 / exp(eta))
d$frac <- rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8)
d$bin <- rbinom(n, 1, plogis(eta))
d$zip <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$hp <- ifelse(d$cnt == 0, 1L, d$cnt)
d$hp[runif(n) < 0.3] <- 0L
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
comp <- rbinom(n, 1, 0.5)
d$mixc <- rpois(n, ifelse(comp == 1, exp(1.5 + 0.7 * d$lat), exp(0.2 + 0.2 * d$lat)))
d$y <- eta + rnorm(n, 0, 0.5)
specs <- list(
  gaussian = c("bf(y ~ me(xo, sx))", "gaussian()"),
  poisson = c("bf(cnt ~ me(xo, sx))", "poisson()"),
  negbin = c("bf(nb ~ me(xo, sx))", "negbinomial()"),
  lognormal = c("bf(pos ~ me(xo, sx))", "lognormal()"),
  gamma = c("bf(gam ~ me(xo, sx))", "Gamma('log')"),
  beta = c("bf(frac ~ me(xo, sx))", "Beta()"),
  bernoulli = c("bf(bin ~ me(xo, sx))", "bernoulli()"),
  zip = c("bf(zip ~ me(xo, sx))", "zero_inflated_poisson()"),
  hurdle = c("bf(hp ~ me(xo, sx))", "hurdle_poisson()"),
  cumulative = c("bf(ord ~ me(xo, sx))", "cumulative()"),
  categorical = c("bf(cat ~ me(xo, sx))", "categorical()"))
nd <- data.frame(xo = c(0, 1, 1, 1), sx = c(0.3, 0.3, 0.8, 1.2))
rowsof <- function(a) if (length(dim(a)) == 3L) a[, , 1] else a
out <- list()
for (nm in names(specs)) {
  cat("==", nm, "\n")
  pr <- tryCatch(rev_pair(specs[[nm]][1], specs[[nm]][2], d, ndraws = 4000), error = function(e) e)
  if (inherits(pr, "error")) {
    cat("  pair ERR", conditionMessage(pr), "\n")
    next
  }
  fit <- pr$fit
  be <- rowsof(posterior_epred(pr$fx, newdata = nd))
  fv <- tryCatch(quiet(fitted(fit, newdata = nd)), error = function(e) e)
  if (inherits(fv, "error")) {
    cat("  fitted ERR", conditionMessage(fv), "\n")
    next
  }
  if (length(dim(fv)) == 3L) fv <- fv[, , 1]
  out[[nm]] <- data.frame(
    family = nm, row = seq_len(nrow(nd)), sx = nd$sx,
    brms_mean = round(colMeans(be), 4), frm_est = round(fv[, "Estimate"], 4),
    ratio = round(colMeans(be) / fv[, "Estimate"], 3),
    brms_sd = round(apply(be, 2, sd), 4), frm_se = round(fv[, "Est.Error"], 4),
    brms_lo = round(apply(be, 2, quantile, 0.025), 4), frm_lo = round(fv[, "Q2.5"], 4),
    brms_hi = round(apply(be, 2, quantile, 0.975), 4), frm_hi = round(fv[, "Q97.5"], 4))
  em <- tryCatch(quiet(as.data.frame(emmeans::emmeans(fit, ~ xo, at = list(xo = 1, sx = 0.8), epred = TRUE))), error = function(e) e)
  emb <- tryCatch(quiet(as.data.frame(emmeans::emmeans(pr$fx, ~ xo, at = list(xo = 1, sx = 0.8), epred = TRUE))), error = function(e) e)
  cat("  emmeans epred frmtmb:", if (inherits(em, "error")) substr(conditionMessage(em), 1, 110) else paste(round(em[[2]], 4), "SE", round(em$SE, 4)), "\n")
  cat("  emmeans epred brms  :", if (inherits(emb, "error")) substr(conditionMessage(emb), 1, 110) else paste(round(emb[[2]], 4), "HPD", round(emb$lower.HPD, 4), round(emb$upper.HPD, 4)), "\n")
  mp <- tryCatch(quiet(marginaleffects::predictions(fit, newdata = nd)), error = function(e) e)
  cat("  predictions se:", if (inherits(mp, "error")) substr(conditionMessage(mp), 1, 110) else round(mp$std.error, 4), "\n")
  ct <- tryCatch(quiet(frm_curve(fit, newdata = nd, transform = TRUE, simultaneous = FALSE)), error = function(e) e)
  cat("  frm_curve transform:", if (inherits(ct, "error")) substr(conditionMessage(ct), 1, 110) else paste(round(ct$.estimate, 3), collapse = " "), "\n")
  ce <- tryCatch(quiet(conditional_effects(fit, "xo", resolution = 2, int_conditions = list(xo = 1), conditions = data.frame(sx = 1.2))[[1]]), error = function(e) e)
  cat("  CE est/lower/upper at sx 1.2:", if (inherits(ce, "error")) substr(conditionMessage(ce), 1, 110) else round(c(ce$estimate__[1], ce$lower__[1], ce$upper__[1]), 3), "\n")
}
res <- do.call(rbind, out)
rownames(res) <- NULL
print(res, row.names = FALSE)
write.table(res, "dev/review-predread-out/mefam.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
