# Review of lane readers: fits with no me(), no latent residual and no new
# level. fitted() in-sample and at new rows of seen levels, predict() (seed
# 1), residuals() and conditional_effects() on both arms, saved to
# dev/review-readers-out/regress-<arm>.rds for dev/review-readers-regress-cmp.R.
# Data seed 9301. Both arms are run pinned to cores 0-1 (start /affinity 3),
# per dev/rtmb-pitfalls.md item 21.
source("dev/review-readers-setup.R")
set.seed(9301)
n <- 400; J <- 15
d <- data.frame(x = rnorm(n), z = runif(n), g = factor(sample(J, n, TRUE)))
u <- rnorm(J, 0, 0.5)
eta <- 0.3 + 0.5 * d$x + u[d$g]
d$gau <- eta + rnorm(n, 0, 0.7)
d$gau2 <- 0.5 * d$gau + rnorm(n)
d$cnt <- rpois(n, exp(eta))
d$nb <- rnbinom(n, mu = exp(eta), size = 2)
d$bin <- rbinom(n, 1, plogis(eta))
d$tr <- 10L
d$bnm <- rbinom(n, 10, plogis(eta))
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                        labels = FALSE))
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$pos <- exp(eta + rnorm(n, 0, exp(-0.5 + 0.3 * d$x)))
lam <- exp(eta)
cc <- rpois(n, lam)
while (any(cc == 0)) cc[cc == 0] <- rpois(sum(cc == 0), lam[cc == 0])
d$hur <- ifelse(runif(n) < plogis(-0.5 + 0.6 * d$x), 0L, cc)
d$zinb <- ifelse(runif(n) < 0.2, 0L, d$nb)
d$gam <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$bet <- plogis(eta + rnorm(n, 0, 0.5))
d$stu <- eta + rt(n, 4) * 0.6
d$wei <- rweibull(n, 2, exp(eta))
d$nly <- 3 * exp(-0.8 * d$z * 3) + rnorm(n, 0, 0.1)
d$smy <- sin(2 * d$x) + rnorm(n, 0, 0.3)
specs <- list(
  gaussian = list(bf(gau ~ x), gaussian()),
  gaussian_re = list(bf(gau ~ x + (1 | g)), gaussian()),
  poisson_re = list(bf(cnt ~ x + (1 | g)), poisson()),
  negbinomial = list(bf(nb ~ x), negbinomial()),
  bernoulli_slope = list(bf(bin ~ x + (1 + x | g)), bernoulli()),
  binomial = list(bf(bnm | trials(tr) ~ x), binomial()),
  cumulative = list(bf(ord ~ x), cumulative()),
  categorical = list(bf(cat ~ x), categorical()),
  lognormal = list(bf(pos ~ x, sigma ~ x), lognormal()),
  hurdle = list(bf(hur ~ x, hu ~ x), hurdle_poisson()),
  zinb = list(bf(zinb ~ x), zero_inflated_negbinomial()),
  gamma = list(bf(gam ~ x), Gamma(link = "log")),
  beta = list(bf(bet ~ x), Beta()),
  student = list(bf(stu ~ x), student()),
  weibull = list(bf(wei ~ x), weibull()),
  mvbf = list(mvbf(bf(gau ~ x), bf(gau2 ~ x)) + set_rescor(TRUE), NULL),
  nl = list(bf(nly ~ a * exp(-b * z), a ~ 1, b ~ 1, nl = TRUE), gaussian()),
  smooth = list(bf(smy ~ s(x)), gaussian()))
nd <- data.frame(x = c(-1, 0, 1), z = c(0.2, 0.5, 0.8), g = d$g[1:3],
                 tr = 10L)
grab <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
out <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  fit <- if (is.null(sp[[2]])) frm(sp[[1]], data = d) else
    frm(sp[[1]], family = sp[[2]], data = d)
  ce_args <- list(fit, effects = if (nm == "nl") "z" else "x")
  if (nm == "nl") ce_args$int_conditions <- list(z = c(0.2, 0.8))
  else ce_args$int_conditions <- list(x = c(-1, 0, 1))
  out[[nm]] <- list(
    fitted_in = grab(fitted(fit)),
    fitted_nd = grab(fitted(fit, newdata = nd)),
    linpred_se = grab(frm_linpred(fit, newdata = nd, se.fit = TRUE)),
    predict = grab({ set.seed(1); predict(fit, newdata = nd) }),
    resid = grab(residuals(fit)),
    ce = grab(lapply(do.call(conditional_effects, ce_args), as.data.frame)))
  cat(nm, "done\n")
}
saveRDS(out, file.path("dev/review-readers-out",
                       paste0("regress-", arm, ".rds")))
