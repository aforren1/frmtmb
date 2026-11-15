# Lane priors2: every bounded fit of a grid, on one build, as one line
# per case (status, code, objective, every outer parameter at 17
# digits), so that two builds can be compared (dev/priors2-survey-cmp.R).
# PRIORS2_LIB unset = lane library; PRIORS2_LIB=none = rellib-r21 only.
# Rscript dev/priors2-survey.R dev/priors2-out/survey-<arm>.tsv
lib <- Sys.getenv("PRIORS2_LIB", "C:/Users/adf44/source/r/rellib-r22")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
out <- commandArgs(trailingOnly = TRUE)[1]
set.seed(1)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:15, 10)),
                t = factor(rep(1:5, 30)))
d$y <- 1 + 0.4 * d$x - 0.3 * d$z + rnorm(15)[d$g] * 0.5 + rnorm(n, sd = 0.7)
d$cnt <- rpois(n, exp(0.3 + 0.3 * d$x))
d$zcnt <- d$cnt * rbinom(n, 1, 0.7)
d$p <- plogis(0.2 + 0.5 * d$x + rnorm(n, sd = 0.5))
d$o <- cut(d$y, 4, labels = FALSE)
d$ys <- 1 + rnorm(n, sd = exp(0.2 + 0.3 * d$x))
d$yt <- 1 + 0.5 * d$x + rt(n, df = 4) * 0.6
sp <- set_prior
# the covariance class is "theta" on 0.82.0
cov_cls <- if ("thetacov" %in% frmtmb:::frmtmb_prior_classes) "thetacov" else
  "theta"
cases <- list(
  # interior maxima under bounds that do not bind
  b_lb_in = list(y ~ x + z, gaussian(), sp("", class = "b", coef = "x",
                                            lb = 0)),
  b_ub_in = list(y ~ x + z, gaussian(), sp("", class = "b", coef = "z",
                                            ub = 0)),
  b_box_in = list(y ~ x + z, gaussian(), sp("normal(0, 1)", class = "b",
                                             lb = -2, ub = 2)),
  b_lognormal_in = list(y ~ x + z, gaussian(),
                        sp("lognormal(0, 1)", class = "b", coef = "x",
                           lb = 0)),
  b_gamma2_in = list(y ~ x + z, gaussian(),
                     sp("gamma(2, 2)", class = "b", coef = "x", lb = 0)),
  b_beta_in = list(y ~ x + z, gaussian(),
                   sp("beta(2, 2)", class = "b", coef = "x", lb = 0,
                      ub = 1)),
  b_invgamma_in = list(y ~ x + z, gaussian(),
                       sp("inv_gamma(3, 1)", class = "b", coef = "x",
                          lb = 0)),
  b_unif_in = list(y ~ x + z, gaussian(),
                   sp("uniform(0, 1)", class = "b", coef = "x")),
  sd_lb_in = list(y ~ x + (1 | g), gaussian(),
                  sp("", class = "sd", lb = 0.1)),
  sd_ub_in = list(y ~ x + (1 | g), gaussian(),
                  sp("", class = "sd", ub = 3)),
  sigma_lb_in = list(y ~ x + z, gaussian(),
                     sp("", class = "sigma", lb = 0.35)),
  sigma_unif_in = list(y ~ x + z, gaussian(),
                       sp("uniform(0.35, 5)", class = "sigma")),
  zi_ub_in = list(zcnt ~ x, zero_inflated_poisson(),
                  sp("", class = "zi", ub = 0.6)),
  zi_unif_in = list(zcnt ~ x, zero_inflated_poisson(),
                    sp("uniform(0, 0.6)", class = "zi")),
  phi_lb_in = list(p ~ x, Beta(), sp("", class = "phi", lb = 0.5)),
  shape_lb_in = list(cnt ~ x, negbinomial(),
                     sp("", class = "shape", lb = 0.3)),
  nu_lb_in = list(yt ~ x, student(), sp("", class = "nu", lb = 1.5)),
  icpt_y1 = list(y ~ 1, gaussian(), sp("", class = "Intercept", lb = 0)),
  nlpar_lb_in = list(bf(y ~ a + b * x, a ~ 1, b ~ 1, nl = TRUE),
                     gaussian(), sp("normal(0, 5)", nlpar = "b", lb = 0) +
                       sp("normal(0, 5)", nlpar = "a")),
  sigma_dpar_b_in = list(bf(ys ~ 1, sigma ~ x), gaussian(),
                         sp("", class = "b", coef = "x", dpar = "sigma",
                            lb = 0)),
  thetacov_in = list(y ~ toep(0 + t | g), gaussian(),
                     sp("", class = cov_cls, coef = "theta_6",
                        lb = -3, ub = 3)),
  # maxima on the bound
  b_lb_on = list(y ~ x + z, gaussian(), sp("", class = "b", coef = "z",
                                            lb = 0)),
  b_ub_on = list(y ~ x + z, gaussian(), sp("", class = "b", coef = "x",
                                            ub = 0.1)),
  b_unif_on = list(y ~ x + z, gaussian(),
                   sp("uniform(0, 0.2)", class = "b", coef = "x")),
  b_exp_on = list(y ~ x + z, gaussian(),
                  sp("exponential(1)", class = "b", coef = "z", lb = 0)),
  sd_lb_on = list(y ~ x + (1 | g), gaussian(),
                  sp("", class = "sd", lb = 1.5)),
  sigma_lb_on = list(y ~ x + z, gaussian(),
                     sp("", class = "sigma", lb = 1.2)),
  sigma_unif_on = list(y ~ x + z, gaussian(),
                       sp("uniform(1.2, 5)", class = "sigma")),
  zi_unif_on = list(zcnt ~ x, zero_inflated_poisson(),
                    sp("uniform(0, 0.1)", class = "zi")),
  phi_ub_on = list(p ~ x, Beta(), sp("", class = "phi", ub = 1)),
  # the data push past a bound where the density is 0
  b_lognormal_push = list(y ~ x + z, gaussian(),
                          sp("lognormal(0, 1)", class = "b", coef = "z",
                             lb = 0)),
  b_gamma2_push = list(y ~ x + z, gaussian(),
                       sp("gamma(2, 2)", class = "b", coef = "z", lb = 0)),
  b_beta_push = list(y ~ x + z, gaussian(),
                     sp("beta(2, 2)", class = "b", coef = "z", lb = 0,
                        ub = 1))
)
fmt <- function(v) paste(sprintf("%.17g", v), collapse = ",")
lines <- character(0)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(suppressWarnings(suppressMessages(
    frm(cs[[1]], data = d, family = cs[[2]], prior = cs[[3]]))),
    error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "error")) {
    lines <- c(lines, paste(nm, "ERROR", NA, NA,
                            gsub("[\t\n]", " ",
                                 substr(conditionMessage(r), 1, 120)),
                            sprintf("%.2f", el), sep = "\t"))
  } else {
    lines <- c(lines, paste(nm, "OK", r$opt$convergence,
                            sprintf("%.17g", r$opt$objective),
                            fmt(r$opt$par), sprintf("%.2f", el), sep = "\t"))
  }
}
writeLines(c(paste("frmtmb", dirname(find.package("frmtmb"))), lines), out)
