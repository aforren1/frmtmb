# Fix round 1, F2: under a likelihood that does not depend on the bounded
# parameter, frm_sample()'s default prior must leave its draws uniform on
# the interval (brms's flat prior on rho). Rscript <this>
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.sample))
lo <- -0.4
hi <- 0.9
frmtmb::frmtmb_register_autocor(
  "flatphi", args = "V",
  build = function(spec, rows, n_all, lookup) list(n = length(rows)),
  npar = 1L, start = function(aux) 0,
  natural = function(theta, aux) c(toyphi = lo + (hi - lo) / (1 + exp(-theta[1]))),
  bounds = function(aux) cbind(lower = lo, upper = hi),
  whiten = function(y, mu, sigma, theta, aux) {
    list(u = (y - mu) / sigma, logjac = -length(y) * log(sigma))
  },
  unwhiten = function(u, mu, sigma, theta, aux) mu + sigma * u,
  sigma_varies = FALSE)
set.seed(3)
d <- data.frame(x = rnorm(30))
d$y <- 1 + d$x + rnorm(30)
t0 <- proc.time()
fit <- withCallingHandlers(
  frm(y ~ x + flatphi(V), data = d, data2 = list(V = 1)),
  warning = function(w) { cat("FRM WARNING:", conditionMessage(w), "\n");
    invokeRestart("muffleWarning") })
print(frmtmb::frmtmb_autocor_pars(fit))
ds <- withCallingHandlers(
  frm_sample(fit, chains = 2, iter = 2000, warmup = 500, refresh = 0,
             seed = 5),
  message = function(m) cat("MSG:", conditionMessage(m)),
  warning = function(w) { cat("SAMPLE WARNING:", conditionMessage(w), "\n");
    invokeRestart("muffleWarning") })
m <- as.matrix(ds)
print(colnames(m))
z <- m[, "toyphi"]
cat("range", range(z), "\n")
cat("mean", mean(z), "want", (lo + hi) / 2, "\n")
cat("sd", sd(z), "want", (hi - lo) / sqrt(12), "\n")
print(quantile(z, c(0.1, 0.25, 0.5, 0.75, 0.9)))
print(lo + (hi - lo) * c(0.1, 0.25, 0.5, 0.75, 0.9))
cat("ess", posterior::ess_basic(z), "\n")
print(ks.test(z, "punif", lo, hi))
cat("time", (proc.time() - t0)[3], "\n")
