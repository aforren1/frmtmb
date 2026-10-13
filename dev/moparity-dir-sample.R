# frm_sample() against brms 2.23.0: a simo prior dirichlet(c(5, 1, 1)),
# a group-level mo() term, and mo() crossed with a two-level factor
# (brms 2.23.0 is right for two levels, brms#1828). 4 chains of 2000
# (1000 warmup), seed 1, on both packages; per weight and coefficient
# the mean, sd, their Monte Carlo errors, R-hat, bulk ESS; divergences.
#   Rscript dev/moparity-dir-sample.R lane|brms
args <- commandArgs(TRUE)
arm <- args[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
out <- "dev/moparity-log"
set.seed(1); n <- 100
x <- sample(0:3, n, TRUE)
weak <- data.frame(x = x, y = 0.05 * x + rnorm(n))
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
strong <- data.frame(income, ls)
set.seed(4); n <- 300; ng <- 20
re <- data.frame(x = sample(0:3, n, TRUE),
                 g = factor(sample(paste0("g", 1:ng), n, TRUE)),
                 f2 = factor(sample(c("a", "b"), n, TRUE)))
bg <- rnorm(ng, 0, 0.5)
re$y <- 1 + (1 + bg[re$g] + 0.8 * (re$f2 == "b")) *
  c(0, 1, 1.5, 1.7)[re$x + 1] + rnorm(n, 0, 0.7)
cases <- list(
  weak = list(f = y ~ mo(x), d = weak, coef = "mox1",
              cols = "^simo_mox1\\["),
  strong = list(f = ls ~ mo(income), d = strong, coef = "moincome1",
                cols = "^simo_moincome1\\["),
  mo_re = list(f = y ~ mo(x) + (mo(x) | g), d = re, coef = NULL,
               cols = "^(simo_mox1\\[|bsp_mox$|sd_g__mox$|sd_g__Intercept$)"),
  mo_fac = list(f = y ~ mo(x) * f2, d = re, coef = NULL,
                cols = "^(simo_|bsp_)"))
summ <- function(M, nch, cols) {
  arr <- array(M[, cols], c(nrow(M) / nch, nch, length(cols)),
               dimnames = list(NULL, NULL, cols))
  s <- posterior::summarise_draws(posterior::as_draws_array(arr),
                                  "mean", "sd", "rhat", "ess_bulk",
                                  "mcse_mean", "mcse_sd")
  as.data.frame(s)
}
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  if (arm == "brms") {
    pr <- if (!is.null(cs$coef)) {
      brms::set_prior("dirichlet(c(5, 1, 1))", class = "simo",
                      coef = cs$coef)
    }
    b <- brms::brm(cs$f, data = cs$d, chains = 4, iter = 2000, seed = 1,
                   cores = 1, refresh = 0, prior = pr)
    M <- as.matrix(b)
    sp <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
  } else {
    pr <- if (!is.null(cs$coef)) {
      set_prior("dirichlet(c(5, 1, 1))", class = "simo", coef = cs$coef)
    }
    fit <- suppressMessages(suppressWarnings(
      frm(bf(cs$f), data = cs$d, family = gaussian())))
    s <- suppressMessages(frm_sample(fit, chains = 4, iter = 2000,
                                     warmup = 1000, seed = 1, cores = 1,
                                     refresh = 0, prior = pr))
    M <- s$draws
    # frm_sample() keeps a block's raw theta in the draws; the us
    # block's first two entries are its log standard deviations
    if (identical(nm, "mo_re")) {
      M <- cbind(M, sd_g__Intercept = exp(M[, "theta_1"]),
                 sd_g__mox = exp(M[, "theta_2"]))
    }
    sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
  }
  cols <- grep(cs$cols, colnames(M), value = TRUE)
  S <- summ(M, 4, cols)
  div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  cat("\n==", nm, arm, ": divergences", div, "\n")
  print(S, digits = 4, row.names = FALSE)
  res[[nm]] <- list(S = S, div = div)
}
saveRDS(res, file.path(out, paste0("dir-sample-", arm, ".rds")))
