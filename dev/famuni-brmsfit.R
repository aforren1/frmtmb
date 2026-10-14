# Fit agreement with brms 2.23.0 on the same simulated data: frmtmb's ML
# estimates against brms's posterior medians, in units of the posterior
# SD, and frmtmb's standard errors against the posterior SDs. n = 500 so
# that brms's default priors carry little weight. One data set per
# family, seed 20261013 + family index. brms: 4 chains x 2000 iterations
# (1000 warmup), seed 1, rstan backend.
# Output: dev/famuni-out/brmsfit.txt; draws in dev/famuni-out/brmsfit-*.rds
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(pkgbuild::has_build_tools())
suppressMessages(library(brms))
n <- 500
sim <- list()
set.seed(20261014)
d <- data.frame(x = rnorm(n))
nu <- 4
d$y <- exp(1 + 0.4 * d$x) / gamma(1 - 1 / nu) * rexp(n)^(-1 / nu)
sim$frechet <- list(d = d, bform = bf(y ~ x), bfam = brms::frechet(),
                    truth = c(Intercept = 1, x = 0.4, nu = nu))
set.seed(20261015)
d <- data.frame(x = rnorm(n), z = rnorm(n))
xi <- 0.2
sg <- exp(-0.2 + 0.3 * d$z)
d$y <- 1 + 0.6 * d$x + sg * (rexp(n)^(-xi) - 1) / xi
sim$gev <- list(d = d, bform = bf(y ~ x, sigma ~ z),
                bfam = brms::gen_extreme_value(),
                truth = c(Intercept = 1, x = 0.6, sigma_Intercept = -0.2,
                          sigma_z = 0.3, xi = xi))
set.seed(20261016)
d <- data.frame(x = rnorm(n), z = rnorm(n))
q <- plogis(1.2 + 0.5 * d$x)
k <- exp(0.3 + 0.25 * d$z)
d$y <- pmax(ceiling((log1p(-runif(n)) / log(q))^(1 / k) - 1), 0)
sim$dw <- list(d = d, bform = bf(y ~ x, shape ~ z),
               bfam = brms::brmsfamily("discrete_weibull"),
               truth = c(Intercept = 1.2, x = 0.5, shape_Intercept = 0.3,
                         shape_z = 0.25))
fits <- list()
for (nm in names(sim)) {
  s <- sim[[nm]]
  t0 <- proc.time()[["elapsed"]]
  fits[[nm]] <- suppressWarnings(brm(s$bform, data = s$d, family = s$bfam,
                                     chains = 4, iter = 2000, cores = 4,
                                     seed = 1, refresh = 0, silent = 2))
  cat(nm, "brms fit in", round(proc.time()[["elapsed"]] - t0), "s\n")
  saveRDS(list(draws = as.data.frame(fits[[nm]]), d = s$d, truth = s$truth,
               rhat = brms::rhat(fits[[nm]]),
               neff = brms::neff_ratio(fits[[nm]])),
          file.path("dev/famuni-out", paste0("brmsfit-", nm, ".rds")))
}

suppressMessages(library(frmtmb, warn.conflicts = FALSE))
cat("frmtmb from", find.package("frmtmb"), "\n")
fams <- list(frechet = frmtmb::frechet(),
             gev = frmtmb::gen_extreme_value(),
             dw = frmtmb::discrete_weibull())
forms <- list(frechet = frmtmb::bf(y ~ x), gev = frmtmb::bf(y ~ x, sigma ~ z),
              dw = frmtmb::bf(y ~ x, shape ~ z))
# brms draw name for each frmtmb coefficient
bname <- c(Intercept = "b_Intercept", x = "b_x",
           sigma_Intercept = "b_sigma_Intercept", sigma_z = "b_sigma_z",
           shape_Intercept = "b_shape_Intercept", shape_z = "b_shape_z",
           nu = "nu", xi = "xi")
for (nm in names(sim)) {
  r <- readRDS(file.path("dev/famuni-out", paste0("brmsfit-", nm, ".rds")))
  fit <- frmtmb::frm(forms[[nm]], family = fams[[nm]], data = r$d)
  cf <- frmtmb::fixef(fit)
  est <- cf[, "Estimate"]
  se <- cf[, "Est.Error"]
  names(est) <- names(se) <- rownames(cf)
  # a dpar without a formula is reported on its own scale; its SE by
  # the delta method through its link
  for (dp in intersect(c("nu", "xi"), names(r$truth))) {
    lp <- frmtmb::frm_linpred(fit, dpar = dp, type = "response",
                              se.fit = TRUE)
    est[dp] <- lp$fit[1]
    se[dp] <- lp$se.fit[1]
  }
  cat("\n===", nm, " (max rhat", format(max(r$rhat, na.rm = TRUE),
                                         digits = 4),
      ", min neff ratio", format(min(r$neff, na.rm = TRUE), digits = 3),
      ")\n")
  cat(sprintf("  %-16s %8s %9s %9s %9s %9s %10s %8s\n", "parameter", "truth",
              "ML", "ML SE", "post med", "post SD", "z (ML-med)", "SE/SD"))
  for (p in names(r$truth)) {
    dr <- r$draws[[bname[[p]]]]
    med <- stats::median(dr)
    sdv <- stats::sd(dr)
    cat(sprintf("  %-16s %8.3f %9.4f %9.4f %9.4f %9.4f %10.3f %8.3f\n", p,
                r$truth[[p]], est[[p]], se[[p]], med, sdv,
                (est[[p]] - med) / sdv, se[[p]] / sdv))
  }
  cat("  frmtmb logLik", format(as.numeric(logLik(fit)), digits = 10),
      "; brms log-lik at its posterior median parameters is not computed\n")
}
