# Fix round: the new rounding-sensitive assertions over 20 ulp-perturbed
# responses, y * (1 + 2e-15 * k), k = 1..20.
# Rscript dev/priortab-ulps2.R > dev/priortab-out/ulps2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-priortab-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- 120
d0 <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d0$y <- 1 + 0.3 * d0$z + rnorm(10)[d0$g] * 0.3 + rnorm(n, sd = 0.3)
q <- function(e) suppressWarnings(suppressMessages(e))
out <- NULL
for (k in 1:20) {
  d <- d0
  d$y <- d$y * (1 + 2e-15 * k)
  d$yneg <- d$y - 0.8 * d$z
  a <- q(frm(y ~ z, data = d, prior = set_prior("uniform(0.6, 1)",
                                               class = "sigma")))
  b <- q(frm(y ~ z, data = d, prior = set_prior("uniform(0.6, 1)",
             class = "sigma", lb = 0.6, ub = 1)))
  ia <- q(frm(bf(y ~ 0 + Intercept + z), data = d,
              prior = set_prior("uniform(0, 0.5)", class = "b",
                                coef = "Intercept")))
  ex <- q(frm(yneg ~ z, data = d, prior = set_prior("exponential(1)",
                                                   class = "b", coef = "z")))
  sp <- function(...) set_prior("beta(2, 2)", class = "b", coef = "z", ...)
  m8 <- q(frm(yneg ~ z, data = d, prior = sp(lb = 1e-8, ub = 1)))
  m10 <- q(frm(yneg ~ z, data = d, prior = sp(lb = 1e-10, ub = 1)))
  at0 <- q(frm(yneg ~ z, data = d, prior = sp(lb = 0, ub = 1)))
  out <- rbind(out, c(
    sigma_identical = identical(a$opt$par, b$opt$par),
    sigma_06 = identical(unname(sigma(a)), 0.6),
    int_05 = identical(unname(fixef(ia)["Intercept", "Estimate"]), 0.5),
    exp_0 = identical(unname(fixef(ex)["z", "Estimate"]), 0),
    beta_rel = abs(m8$opt$objective - m10$opt$objective) /
      abs(m8$opt$objective),
    beta_z = unname(fixef(m8)["z", "Estimate"]),
    at0_code = at0$opt$convergence))
}
print(out)
cat("all identical cols TRUE:", all(out[, 1:4] == 1),
    " max beta_rel", max(out[, "beta_rel"]), " min beta_z",
    min(out[, "beta_z"]), " at0 codes", unique(out[, "at0_code"]), "\n")
