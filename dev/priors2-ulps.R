# Lane priors2: the new fit assertions over 20 perturbations of the
# response by a few ulps (y * (1 + 2e-15 * k), k = 0..19), so that an
# assertion is not one rounding's luck.
# Rscript dev/priors2-ulps.R > dev/priors2-out/ulps.txt
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
set.seed(1)
dz <- data.frame(z = rnorm(n))
invisible(rnorm(10)); invisible(rnorm(n))
dz$cnt <- rpois(n, 2) * rbinom(n, 1, 0.5)
set.seed(1)
ds <- data.frame(z = rnorm(n))
ds$y <- 1 + 0.3 * ds$z + rnorm(n, sd = 0.3)
set.seed(3)
dm <- data.frame(y = c(rnorm(48, 0, 1), rnorm(112, 4, 1)))
tally <- list()
add <- function(nm, ok) tally[[nm]] <<- c(tally[[nm]], isTRUE(ok))
for (k in 0:19) {
  m <- 1 + 2e-15 * k
  dk <- transform(d, yneg = yneg * m, y = y * m)
  for (dn in c("lognormal(0, 1)", "beta(2, 2)", "gamma(2, 2)")) {
    ub <- if (startsWith(dn, "beta")) 1 else NA
    f <- q(frm(yneg ~ z, data = dk,
               prior = set_prior(dn, class = "b", coef = "z", lb = 0,
                                 ub = ub)))
    z <- unname(fixef(f)["z", "Estimate"])
    add(paste("open", dn), f$opt$convergence == 0 &&
          is.finite(f$opt$objective) && z > 0)
  }
  f <- q(frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
             prior = set_prior("uniform(0, 0.3)", class = "zi")))
  zi <- plogis(unname(f$estimates$betad[1]))
  add("zi uniform(0, 0.3)", f$opt$convergence == 0 && zi <= 0.3 &&
        (0.3 - zi) / 0.3 < 1e-8)
  f <- q(frm(y ~ z, data = transform(ds, y = y * m),
             prior = set_prior("uniform(0.35, 5)", class = "sigma",
                               lb = 0.35, ub = 5)))
  s <- exp(unname(f$estimates$betad[1]))
  add("sigma uniform(0.35, 5)", f$opt$convergence == 0 && s >= 0.35 &&
        (s - 0.35) / 0.35 < 1e-8)
  f <- q(frm(y ~ 1, data = transform(dm, y = y * m),
             family = mixture(gaussian(), gaussian()),
             prior = set_prior("dirichlet(c(3, 2))", class = "theta")))
  w <- plogis(unname(f$estimates$betad[3]))
  lp <- lgamma(5) - lgamma(3) - lgamma(2) + 3 * log(w) + 2 * log1p(-w)
  add("theta MAP objective = nll - log prior",
      f$opt$convergence == 0 &&
        abs(f$opt$objective + as.numeric(logLik(f)) + lp) <
        1e-8 * abs(f$opt$objective))
}
for (nm in names(tally)) {
  cat(sprintf("%-40s %d of %d\n", nm, sum(tally[[nm]]),
              length(tally[[nm]])))
}
