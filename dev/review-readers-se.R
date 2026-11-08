# Review of lane readers: conditional_effects()'s se__ and band columns on
# both arms, saved to dev/review-readers-out/se-<arm>.rds and compared by
# dev/review-readers-se-cmp.R. Data seed 7971 as the lane's test, plus a
# hurdle (seed 7972), a bernoulli and a cumulative (seed 7973), a gamma with
# an inverse link (seed 7974).
source("dev/review-readers-setup.R")
set.seed(7971)
n <- 300
d <- data.frame(x = rnorm(n))
eta <- 0.4 + 0.7 * d$x
d$cnt <- rpois(n, exp(eta))
d$pos <- exp(eta + rnorm(n, 0, exp(-0.5 + 0.3 * d$x)))
set.seed(7972)
lam <- exp(0.5 + 0.5 * d$x)
cc <- rpois(n, lam)
while (any(cc == 0)) cc[cc == 0] <- rpois(sum(cc == 0), lam[cc == 0])
d$hur <- ifelse(runif(n) < plogis(-0.5 + 0.6 * d$x), 0L, cc)
set.seed(7973)
d$bin <- rbinom(n, 1, plogis(0.3 + 0.9 * d$x))
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                        labels = FALSE))
set.seed(7974)
d$gam <- rgamma(n, shape = 3, rate = 3 * (0.5 + 0.2 * (d$x + 3)))
d$gau <- 1 + 0.5 * d$x + rnorm(n)
ic <- list(x = c(-1, 0, 1))
specs <- list(
  gaussian = list(bf(gau ~ x), gaussian(), NULL, list()),
  poisson = list(bf(cnt ~ x), poisson(), NULL, list()),
  poisson_predict = list(bf(cnt ~ x), poisson(), NULL,
                         list(method = "predict")),
  lognormal = list(bf(pos ~ x, sigma ~ x), lognormal(), NULL, list()),
  lognormal_sigma = list(bf(pos ~ x, sigma ~ x), lognormal(), "sigma",
                         list()),
  lognormal_pp = list(bf(pos ~ x, sigma ~ x), lognormal(), NULL,
                      list(method = "posterior_predict")),
  hurdle = list(bf(hur ~ x, hu ~ x), hurdle_poisson(), NULL, list()),
  hurdle_hu = list(bf(hur ~ x, hu ~ x), hurdle_poisson(), "hu", list()),
  bernoulli = list(bf(bin ~ x), bernoulli(), NULL, list()),
  cumulative_cat = list(bf(ord ~ x), cumulative(), NULL,
                        list(categorical = TRUE)),
  cumulative_mean = list(bf(ord ~ x), cumulative(), NULL, list()),
  gamma_inverse = list(bf(gam ~ x), Gamma(link = "inverse"), NULL, list()),
  poisson_profile = list(bf(cnt ~ x), poisson(), NULL,
                         list(band = "profile")))
out <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  fit <- frm(sp[[1]], family = sp[[2]], data = d)
  r <- tryCatch({
    ce <- do.call(conditional_effects,
                  c(list(fit, effects = "x", dpar = sp[[3]],
                         int_conditions = ic), sp[[4]]))[[1]]
    keep <- intersect(c("x", "cats__", "estimate__", "se__", "lower__",
                        "upper__"), names(ce))
    ce[, keep]
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  out[[nm]] <- list(ce = r,
                    fitted = tryCatch(fitted(fit, newdata = data.frame(
                      x = ic$x), dpar = sp[[3]]),
                      error = function(e) NULL))
  cat("\n==", nm, "\n")
  print(r)
}
saveRDS(out, file.path("dev/review-readers-out",
                       paste0("se-", arm, ".rds")))
