# Review of lane loglik, attack 4: the nonlinear restart messages under a
# prior. Does the "log-likelihood" they print match logLik() of the
# returned fit, or logLik_penalized()? Seeds 501 to 503.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
grab <- function(expr) {
  msgs <- character(0)
  f <- withCallingHandlers(expr,
    message = function(m) { msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage") },
    warning = function(w) { msgs <<- c(msgs, paste("W:", conditionMessage(w)))
      invokeRestart("muffleWarning") })
  list(f = f, msgs = msgs)
}
for (seed in 501:503) {
  set.seed(seed)
  d <- data.frame(x = runif(80, 0, 3))
  d$y <- 2 * exp(-0.8 * d$x) + rnorm(80, 0, 0.2)
  # a and b start at the prior locations, 0, so the body a * b * x is 0 on
  # every row and both have no gradient there: a saddle, the restart path
  d$y <- 0.8 * d$x + rnorm(80, 0, 0.3)
  r <- grab(frm(bf(y ~ a * b * x, a + b ~ 1, nl = TRUE), data = d,
                prior = set_prior("normal(0, 1)", nlpar = "a") +
                  set_prior("normal(0, 1)", nlpar = "b")))
  hit <- grep("log-likelihood", r$msgs, value = TRUE)
  cat(sprintf("seed %d  logLik %.6f  logLik_penalized %.6f\n", seed,
              as.numeric(logLik(r$f)), as.numeric(logLik_penalized(r$f))))
  for (m in hit) {
    nums <- regmatches(m, gregexpr("log-likelihood -?[0-9.]+", m))[[1]]
    cat("   message numbers:", paste(nums, collapse = " | "), "\n")
  }
  if (!length(hit)) cat("   no restart message\n")
}
cat("done\n")
