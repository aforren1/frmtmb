# Reviewer: the same newdata check through frmtmb.sample's draws-level
# posterior_predict() (binomial trials()). Seed 7.
# Log dev/bin2port-rev-out/newdata2.txt. FRMTMB_STAN_CACHE must be set.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(7)
n <- 200
db <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
db$y <- rbinom(n, db$N, plogis(0.3 + 0.5 * db$x))
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = db)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fb, chains = 1, iter = 300, warmup = 150, refresh = 0,
             seed = 7)))
ndb <- data.frame(x = c(0, 0), N = c(1000, 2))
p <- tryCatch(posterior_predict(ds, newdata = ndb),
              error = function(e) conditionMessage(e))
if (is.matrix(p)) {
  cat("posterior_predict(draws, newdata N = 1000, 2): column means",
      format(colMeans(p), digits = 4), " max", apply(p, 2, max), "\n")
} else cat("posterior_predict(draws, newdata):", p, "\n")
e <- posterior_epred(ds, newdata = ndb)
cat("posterior_epred(draws, newdata): column means",
    format(colMeans(e), digits = 4), "\n")
p0 <- posterior_predict(ds)
cat("in-sample: mean of draws / N, first 3 rows:",
    format(colMeans(p0)[1:3] / db$N[1:3], digits = 3), "\n")
