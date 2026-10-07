# Lane optima, punch round 1, B1: frm_sample() on a mo() model against
# brms 2.23.0 with its default dirichlet(1). The reviewer's two cases
# (dev/optima-rev-sample2-arm.R): a weakly identified simplex
# (y = 0.05 x + noise, n = 100, D = 3) and brms_monotonic's fit1.
# 4 chains, iter 2000, warmup 1000, seed 1, on both packages. Per
# simplex weight: mean and sd on each side, the difference over its
# Monte Carlo standard error, R-hat and bulk ESS; divergences.
#   Rscript dev/optima-p1-sample.R lane|base [brms]
args <- commandArgs(TRUE)
arm <- args[1]
do_brms <- length(args) > 1 && args[2] == "brms"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "; frmtmb.sample from",
    find.package("frmtmb.sample"), "\n")
out <- "dev/optima-log"
set.seed(1); n <- 100
x <- sample(0:3, n, TRUE)
weak <- data.frame(x = x, y = 0.05 * x + rnorm(n))
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
strong <- data.frame(income, ls)
cases <- list(weak = list(f = y ~ mo(x), d = weak, w = "simo_mox1"),
              strong = list(f = ls ~ mo(income), d = strong,
                            w = "simo_moincome1"))
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
  if (do_brms) {
    b <- brms::brm(cs$f, data = cs$d, chains = 4, iter = 2000, seed = 1,
                   cores = 1, refresh = 0)
    M <- as.matrix(b)
    sp <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
    nch <- 4
  } else {
    fit <- frm(bf(cs$f), data = cs$d, family = gaussian())
    s <- suppressMessages(frm_sample(fit, chains = 4, iter = 2000,
                                     warmup = 1000, seed = 1, cores = 1,
                                     refresh = 0))
    M <- s$draws
    sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
    nch <- 4
    cat("\n==", nm, "draws columns:", paste(colnames(M), collapse = " "),
        "\n")
  }
  cols <- grep(paste0("^", cs$w, "\\["), colnames(M), value = TRUE)
  S <- summ(M, nch, cols)
  div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  cat("\n==", nm, if (do_brms) "brms" else arm, ": divergences", div, "\n")
  print(S, digits = 4, row.names = FALSE)
  res[[nm]] <- list(S = S, div = div)
}
saveRDS(res, file.path(out, paste0("p1-sample-",
                                   if (do_brms) "brms" else arm, ".rds")))
