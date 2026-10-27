# Reviewer: where are frm_sample()'s divergent transitions on the cosy
# data of dev/acprior-agree.R (seed 3, 4 chains of 2000, random init)?
# Log: dev/review-acprior-log/cosydiv2.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(20261027)
x0 <- rnorm(48)
for (i in 1:48) invisible(rnorm(1))
ng2 <- 30; k2 <- 4
dc <- data.frame(t = rep(seq_len(k2), ng2), g = rep(seq_len(ng2), each = k2),
                 x = rnorm(ng2 * k2))
dc$y <- 1 + 0.5 * dc$x + rnorm(ng2, 0, sqrt(0.03))[dc$g] +
  rnorm(nrow(dc), 0, sqrt(0.97))
fit <- frm(y ~ x + cosy(time = t, gr = g), data = dc)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 4, iter = 2000, seed = 3, refresh = 0,
             init = "random")))
cat("class:", class(ds), "\n")
sf <- ds$stanfit %||% ds$fit %||% attr(ds, "stanfit")
if (is.null(sf)) {
  cat("names:", names(ds), "\n")
} else {
  np <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  div <- unlist(lapply(np, function(m) m[, "divergent__"])) == 1
  cosy <- as.matrix(ds)[, "cosy"]
  cat("divergent:", sum(div), "\n")
  cat("cosy at divergent draws:", signif(sort(cosy[div]), 3), "\n")
  cat("share of all draws with cosy below the largest divergent one:",
      mean(cosy <= max(cosy[div])), "\n")
  cat("treedepth max:", max(unlist(lapply(np, function(m) m[, "treedepth__"]))),
      "; stepsize:", signif(sapply(np, function(m) m[1, "stepsize__"]), 3), "\n")
}
