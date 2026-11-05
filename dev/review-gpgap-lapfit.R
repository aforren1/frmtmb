# Short check of fix round 1, item 1: the route the warning leads to. frm() warns,
# so the user calls frm_sample(fit, laplace = TRUE), from frm()'s collapsed sigma
# (the mode start). Does it reach the target? The lane's Laplace runs are the
# formula route (init "random"). And the message's cost at 450 positions. Seed 6601.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-lapfit.R lapfit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
for (n in c(15L, 40L)) {
  set.seed(6601)
  x <- round(runif(n, 0, 6), 3)
  e0 <- rnorm(n)
  d <- data.frame(x = x, y = sin(2 * x) + 1e-3 * e0)
  msg <- NULL
  f0 <- withCallingHandlers(frm(bf(y ~ gp(x)), family = gaussian(), data = d),
                            frmtmb_gp_nugget = function(w) msg <<- conditionMessage(w),
                            warning = function(w) invokeRestart("muffleWarning"))
  cat(sprintf("n %d frm sigma %.3g, warned %s\n", n, as.numeric(sigma(f0)), !is.null(msg)))
  if (n == 15L) cat("MESSAGE:", msg, "\n")
  for (ini in list(NULL, "random")) {
    t0 <- proc.time()[["elapsed"]]
    told <- character()
    args <- list(f0, laplace = TRUE, chains = 4, iter = 2000, warmup = 1000, seed = 21,
                 refresh = 0, cores = 1, control = list(adapt_delta = 0.95))
    if (!is.null(ini)) args$init <- ini
    fl <- withCallingHandlers(tryCatch(do.call(frm_sample, args), error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      NULL
    }), warning = function(w) {
      told <<- c(told, gsub("[\r\n]+", " ", substr(conditionMessage(w), 1, 100)))
      invokeRestart("muffleWarning")
    })
    if (is.null(fl)) next
    dd <- posterior::as_draws_df(fl)
    m <- posterior::extract_variable_matrix(dd, "sigma")
    nm <- grep("^sdgp", posterior::variables(dd), value = TRUE)[1]
    m2 <- posterior::extract_variable_matrix(dd, nm)
    q <- quantile(as.numeric(m), c(.05, .5, .95))
    cat(sprintf(paste0("LAPFIT n %d init %-14s sigma %.3e %.3e %.3e rhat %.3f | sdgp q50 %.3f rhat %.3f | ",
                       "%.0f s\n"), n, if (is.null(ini)) "mode (default)" else ini, q[1], q[2], q[3],
                posterior::rhat(m), median(as.numeric(m2)), posterior::rhat(m2),
                proc.time()[["elapsed"]] - t0))
    cat("  told:", paste(unique(told), collapse = " | "), "\n")
  }
}
# cost at 450 positions
set.seed(6602)
n <- 450
d <- data.frame(x = round(runif(n, 0, 6), 4))
d$y <- sin(2 * d$x) + 0.1 * rnorm(n)
t1 <- proc.time()[["elapsed"]]
f <- suppressWarnings(frm(bf(y ~ gp(x)), data = d))
t2 <- proc.time()[["elapsed"]]
m <- frmtmb:::gp_nugget_message(f)
t3 <- proc.time()[["elapsed"]]
cat(sprintf("COST n %d fit %.1f s, message %.2f s (%.2f percent)\n", n, t2 - t1, t3 - t2,
            100 * (t3 - t2) / (t2 - t1)))
