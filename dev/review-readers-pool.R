# Review of lane readers: conditional_effects() of a frm_multiple() result
# pools by Rubin's rules on ce_pool_scale()'s scale, which the lane changed
# (se_on = "response", the band's own link). Both arms save
# dev/review-readers-out/pool-<arm>.rds; the after arm compares. Data seed
# 7975, three imputations made by jittering x.
source("dev/review-readers-setup.R")
set.seed(7975)
n <- 300
base <- data.frame(x = rnorm(n))
eta <- 0.4 + 0.7 * base$x
base$cnt <- rpois(n, exp(eta))
base$pos <- exp(eta + rnorm(n, 0, exp(-0.5 + 0.3 * base$x)))
lam <- exp(0.5 + 0.5 * base$x)
cc <- rpois(n, lam)
while (any(cc == 0)) cc[cc == 0] <- rpois(sum(cc == 0), lam[cc == 0])
base$hur <- ifelse(runif(n) < plogis(-0.5 + 0.6 * base$x), 0L, cc)
base$bin <- rbinom(n, 1, plogis(eta))
imps <- lapply(1:3, function(i) transform(base, x = x + rnorm(n, 0, 0.1)))
specs <- list(
  poisson = list(bf(cnt ~ x), poisson(), NULL),
  lognormal = list(bf(pos ~ x, sigma ~ x), lognormal(), NULL),
  lognormal_sigma = list(bf(pos ~ x, sigma ~ x), lognormal(), "sigma"),
  hurdle = list(bf(hur ~ x, hu ~ x), hurdle_poisson(), NULL),
  bernoulli = list(bf(bin ~ x), bernoulli(), NULL))
out <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  fm <- frm_multiple(sp[[1]], family = sp[[2]], data = imps)
  out[[nm]] <- tryCatch(
    conditional_effects(fm, effects = "x", dpar = sp[[3]])[[1]][,
      c("x", "estimate__", "se__", "lower__", "upper__")],
    error = function(e) paste("ERROR:", conditionMessage(e)))
}
saveRDS(out, file.path("dev/review-readers-out", paste0("pool-", arm,
                                                        ".rds")))
if (arm == "after" && file.exists("dev/review-readers-out/pool-before.rds")) {
  b <- readRDS("dev/review-readers-out/pool-before.rds")
  for (nm in names(out)) {
    x <- b[[nm]]; y <- out[[nm]]
    if (is.character(x) || is.character(y)) {
      cat(nm, ":", if (is.character(x)) x else "ok", "|",
          if (is.character(y)) y else "ok", "\n")
      next
    }
    cat(sprintf(paste("%-16s pooled est rel %.1e lower rel %.1e upper rel",
                      "%.1e se__ after/before %s\n"), nm,
                max(abs(y$estimate__ / x$estimate__ - 1)),
                max(abs(y$lower__ / x$lower__ - 1)),
                max(abs(y$upper__ / x$upper__ - 1)),
                paste(format(y$se__ / x$se__, digits = 4), collapse = " ")))
  }
}
