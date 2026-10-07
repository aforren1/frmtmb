# Reviewer of lane optima, task 2: frm_sample() posterior of a mo()
# simplex on base vs lane, weakly and strongly identified data.
#   Rscript dev/optima-rev-sample2-arm.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev-out"
ms <- if (exists("mo_simplex", asNamespace("frmtmb")) &&
          arm == "lane") {
  utils::getFromNamespace("mo_simplex", "frmtmb")
} else {
  function(z) { x <- exp(c(0, z)); x / sum(x) }
}
cat("mo_simplex(c(0,0)) =", ms(c(0, 0)), " mo_simplex(c(3,-1)) =",
    format(ms(c(3, -1)), digits = 4), "\n")

set.seed(1); n <- 100
x <- sample(0:3, n, TRUE)
weak <- data.frame(x = x, y = 0.05 * x + rnorm(n))
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
strong <- data.frame(income, ls)
saveRDS(list(weak = weak, strong = strong), file.path(out, "data.rds"))

cases <- list(weak = list(f = bf(y ~ mo(x)), d = weak),
              strong = list(f = bf(ls ~ mo(income)), d = strong))
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(cs$f, data = cs$d, family = gaussian())
  zi <- grep("^zeta", names(fit$estimates), value = TRUE)
  cat("\n==", nm, "== ML zeta", format(fit$estimates[[zi]], digits = 5),
      " ML simplex", format(ms(fit$estimates[[zi]]), digits = 4),
      " logLik", format(as.numeric(logLik(fit)), digits = 8), "\n")
  t0 <- Sys.time()
  s <- withCallingHandlers(
    frm_sample(fit, chains = 2, iter = 2000, warmup = 1000, seed = 1,
               cores = 1, refresh = 0),
    message = function(m) { cat("MSG:", conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { cat("WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  cat("sampling time", format(Sys.time() - t0), "\n")
  sf <- s$stanfit
  a <- rstan::extract(sf, permuted = FALSE)
  cat("stan par names:", paste(dimnames(a)[[3]], collapse = " "), "\n")
  cat("draws colnames:", paste(colnames(s$draws), collapse = " "), "\n")
  am <- as.matrix(s)
  cat("as.matrix colnames:", paste(colnames(am), collapse = " "), "\n")
  sm <- rstan::summary(sf)$summary
  print(round(sm[, c("mean", "sd", "n_eff", "Rhat")], 4))
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  div <- sapply(sp, function(x) sum(x[, "divergent__"]))
  td <- sapply(sp, function(x) max(x[, "treedepth__"]))
  cat("divergent per chain:", div, " max treedepth per chain:", td, "\n")
  res[[nm]] <- list(arr = a, summary = sm, div = div, treedepth = td,
                    ml_zeta = fit$estimates[[zi]], asmat = am,
                    draws = s$draws)
}
saveRDS(res, file.path(out, paste0("draws-", arm, ".rds")))

# does a simo prior row get through?
cat("\n== simo prior row ==\n")
pr <- brms::prior(dirichlet(2), class = simo, coef = moincome1)
r1 <- tryCatch({
  frm(bf(ls ~ mo(income)), data = strong, family = gaussian(), prior = pr)
  "accepted"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("frm(prior = simo row):", r1, "\n")
fit <- frm(bf(ls ~ mo(income)), data = strong, family = gaussian())
r2 <- tryCatch({
  suppressMessages(frm_sample(fit, chains = 1, iter = 50, warmup = 25,
                              seed = 1, cores = 1, refresh = 0, prior = pr))
  "accepted"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("frm_sample(prior = simo row):", r2, "\n")
pr0 <- brms::prior_string("", class = "simo", coef = "moincome1")
r3 <- tryCatch({
  frm(bf(ls ~ mo(income)), data = strong, family = gaussian(), prior = pr0)
  "accepted"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("frm(prior = flat simo row):", r3, "\n")
