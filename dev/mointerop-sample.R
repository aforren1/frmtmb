## frm_sample() of mo() with two values and with one value, on the data
## of dev/mointerop-brms-mo12.R (seed 5, 200 rows), against brms 2.23.0's
## own program on the same data (2 chains of 2000, seed 3 each).
##   Rscript dev/mointerop-sample.R <lane|base> > dev/mointerop-out/sample-<arm>.txt
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "| frmtmb.sample",
    dirname(find.package("frmtmb.sample")), "\n")
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
d$x1 <- 4
w <- character(0)
ds <- withCallingHandlers(tryCatch(suppressMessages(frm_sample(
  y ~ mo(x) + z, data = d, chains = 2, iter = 2000, seed = 3,
  refresh = 0)), error = function(e) e),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
for (m in w) cat("WARNING:", substr(m, 1, 200), "\n")
if (inherits(ds, "error")) {
  cat("frm_sample ERROR:", conditionMessage(ds), "\n")
} else {
  cat("variables:", variables(ds), "\n")
  mf <- as.matrix(ds)
  cat("simo_mox1[1] range:", range(mf[, "simo_mox1[1]"]), "\n")
  code <- brms::stancode(y ~ mo(x) + z, data = d)
  sdat <- brms::standata(y ~ mo(x) + z, data = d)
  mod <- rstan::stan_model(model_code = code)
  sf <- rstan::sampling(mod, data = sdat, chains = 2, iter = 2000,
                        seed = 3, refresh = 0)
  mb <- as.matrix(sf, pars = c("b_Intercept", "b", "bsp", "sigma"))
  colnames(mb) <- c("b_Intercept", "b_z", "bsp_mox", "sigma")
  mc <- function(v) posterior::mcse_mean(matrix(v, ncol = 2))
  cat(sprintf("%-12s %10s %10s %10s %10s %7s %7s\n", "par", "frm mean",
              "frm sd", "brms mean", "brms sd", "z", "sd rat"))
  for (p in colnames(mb)) {
    z <- (mean(mf[, p]) - mean(mb[, p])) /
      sqrt(mc(mf[, p])^2 + mc(mb[, p])^2)
    cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.4f %7.2f %7.3f\n", p,
                mean(mf[, p]), sd(mf[, p]), mean(mb[, p]), sd(mb[, p]), z,
                sd(mf[, p]) / sd(mb[, p])))
  }
  cat("brms max rhat:", max(apply(mb, 2, function(v)
    posterior::rhat(matrix(v, ncol = 2)))), "\n")
  e <- posterior_epred(ds, newdata = data.frame(x = 1:2, z = 0))
  cat("epred(x = 2) - epred(x = 1) mean:", mean(e[, 2] - e[, 1]),
      " bsp_mox mean:", mean(mf[, "bsp_mox"]), "\n")
  r <- tryCatch(posterior_epred(ds, newdata = data.frame(x = 3, z = 0)),
                error = function(e) conditionMessage(e))
  cat("epred x = 3:", if (is.character(r)) r else "returned", "\n")
  print(summary(ds))
  ll <- log_lik(ds)
  cat("log_lik dim:", dim(ll), "finite:", all(is.finite(ll)), "\n")
}
r1 <- tryCatch(suppressMessages(frm_sample(
  y ~ mo(x1) + z, data = d, chains = 1, iter = 200, seed = 3,
  refresh = 0)), error = function(e) e)
cat("\none value:", if (inherits(r1, "error")) {
  paste("ERROR:", conditionMessage(r1))
} else "sampled", "\n")
