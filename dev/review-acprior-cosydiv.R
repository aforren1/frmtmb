# Reviewer: on the cosy data of dev/acprior-agree.R (ML cosy -0.0089),
# frm_sample() after the lane reported 27 divergent transitions and
# brms none. Is that the start at the clamped mode or the bounded
# posterior itself? Mode starts against random starts, three seeds.
# Usage: Rscript dev/review-acprior-cosydiv.R
# Log: dev/review-acprior-log/cosydiv.txt (round 0), r1-cosydiv.txt (fix round 1)
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
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(20261027)
invisible(rnorm(48)); for (i in 1:48) invisible(rnorm(1))
# rebuild exactly as dev/acprior-agree.R's make_data()
make_data <- function() {
  set.seed(20261027)
  ng <- 8; k <- 6
  da <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                   x = rnorm(ng * k))
  e <- numeric(nrow(da))
  for (i in seq_len(nrow(da))) {
    e[i] <- if (da$t[i] == 1) rnorm(1) else 0.9 * e[i - 1] +
      rnorm(1, 0, sqrt(1 - 0.81))
  }
  da$y <- 1 + 0.5 * da$x + e
  ng2 <- 30; k2 <- 4
  dc <- data.frame(t = rep(seq_len(k2), ng2),
                   g = rep(seq_len(ng2), each = k2), x = rnorm(ng2 * k2))
  dc$y <- 1 + 0.5 * dc$x + rnorm(ng2, 0, sqrt(0.03))[dc$g] +
    rnorm(nrow(dc), 0, sqrt(0.97))
  dc
}
dc <- make_data()
fit <- frm(y ~ x + cosy(time = t, gr = g), data = dc)
cat("ML cosy", format(unlist(frmtmb:::autocor_natural(fit$estimates$thetaac,
                                                     fit$frame$autocor[[1]]))),
    "\n")
run <- function(init, seed) {
  w <- character(0)
  ds <- withCallingHandlers(
    suppressMessages(frm_sample(fit, chains = 4, iter = 2000, seed = seed,
                                refresh = 0, init = init)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  dv <- regmatches(w, regexpr("There were [0-9]+ divergent", w))
  cosy <- as.matrix(ds)[, "cosy"]
  cat(sprintf("init %-13s seed %d: %s; clamp warning %s; cosy mean %.4f\n",
              init, seed, if (length(dv)) dv else "0 divergent",
              any(grepl("violates the requested", w)), mean(cosy)))
}
for (s in 3:5) {
  run("last.par.best", s)
  run("random", s)
}
# brms 2.23.0 on the same data and seeds, its own defaults
suppressPackageStartupMessages(library(brms))
fo <- y ~ x + cosy(time = t, gr = g)
code <- stancode(fo, data = dc)
mod <- rstan::stan_model(model_code = code)
for (s in 3:5) {
  sf <- suppressWarnings(rstan::sampling(mod, data = standata(fo, data = dc),
                                         chains = 4, iter = 2000, seed = s,
                                         refresh = 0))
  cat(sprintf("brms seed %d: %d divergent; cosy mean %.4f\n", s,
              rstan::get_num_divergent(sf),
              mean(as.matrix(sf, pars = "cosy"))))
}
