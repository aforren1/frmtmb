# Reviewer: dev/moparity-cox-sample.R at 4 chains of 6000 (1000 warmup), sampler seed 3, data seed 5 only.
# sbhaz[k] and b_x. Seeds 2 and 5 of dev/optima-rev-cox1.R's data (the
# ones whose weights the softmax left at 1e-8), 4 chains of 2000 (1000
# warmup), sampler seed 1, both packages under their default priors.
#   Rscript dev/moparity-cox-sample.R base|lane|brms
args <- commandArgs(TRUE)
arm <- args[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
summ <- function(M, nch, cols) {
  arr <- array(M[, cols], c(nrow(M) / nch, nch, length(cols)),
               dimnames = list(NULL, NULL, cols))
  s <- posterior::summarise_draws(posterior::as_draws_array(arr),
                                  "mean", "sd", "rhat", "ess_bulk",
                                  "mcse_mean", "mcse_sd")
  as.data.frame(s)
}
res <- list()
for (seed in c(5)) {
  set.seed(seed)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
  if (arm == "brms") {
    b <- brms::brm(time | cens(cens) ~ x, data = dd, family = brms::cox(),
                   chains = 4, iter = 6000, warmup = 1000, seed = 3, cores = 1,
                   refresh = 0)
    M <- as.matrix(b)
    sp <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
  } else {
    fit <- suppressWarnings(frm(bf(time | cens(cens) ~ x), family = cox(),
                                data = dd))
    s <- suppressWarnings(suppressMessages(frm_sample(
      fit, chains = 4, iter = 6000, warmup = 1000, seed = 3, cores = 1,
      refresh = 0)))
    M <- s$draws
    sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
    cat("draws columns:", paste(colnames(M), collapse = " "), "\n")
  }
  cols <- grep("^(sbhaz|b_x$)", colnames(M), value = TRUE)
  S <- if (length(cols)) summ(M, 4, cols)
  div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  cat("\n== seed", seed, arm, ": divergences", div, "\n")
  print(S, digits = 4, row.names = FALSE)
  res[[as.character(seed)]] <- list(S = S, div = div)
}
saveRDS(res, file.path("dev/moparity-rev-log", paste0("cox-long-", arm,
                                                  ".rds")))
