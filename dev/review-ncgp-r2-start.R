# Fix round 1, item 5: frm_sample(fit)'s start rule on fresh seeds. Args: seed
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
suppressPackageStartupMessages(library(frmtmb.sample))
seed <- as.integer(args[1])
for (s in c(1e-6, 1e-4, 1e-3)) {
  set.seed(seed); n <- 15L; x <- sort(round(runif(n, 0, 6), 3)); y <- sin(2 * x) + s * rnorm(n)
  fit <- suppressWarnings(frm(y ~ gp(x), data.frame(x = x, y = y)))
  msgs <- character()
  t0 <- proc.time()[[3]]
  dr <- withCallingHandlers(suppressWarnings(frm_sample(fit, chains = 2, iter = 1000, seed = seed, refresh = 0, control = list(max_treedepth = 10))),
                            message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  a <- array(dr$draws[, c("sigma", "sdgp_gpx")], c(nrow(dr$draws) / 2, 2, 2))
  cat(sprintf("RST2 seed %d s %.0e frm sigma %.3g nuggetmsg %s | draws sigma q50 %.3g (%.3g of truth) rhat %.2f sdgp rhat %.2f wall %.0f s\n", seed, s, sigma(fit)[1], !is.null(gp_nugget_message(fit)), median(a[, , 1]), median(a[, , 1]) / s,
              posterior::rhat(a[, , 1]), posterior::rhat(a[, , 2]), proc.time()[[3]] - t0))
  if (any(grepl("chains start at random", msgs))) cat("RST2 MSG", msgs[grepl("chains start", msgs)], "\n")
  flush.console()
  w <- tryCatch(suppressWarnings(suppressMessages(frm_sample(fit, laplace = TRUE, chains = 2, iter = 600, seed = seed, refresh = 0))), error = function(e) e)
  if (inherits(w, "error")) cat("RST2 laplace ERR", conditionMessage(w), "\n") else {
    b <- array(w$draws[, c("sigma", "sdgp_gpx")], c(nrow(w$draws) / 2, 2, 2))
    cat(sprintf("RST2 laplace s %.0e sigma q50 %.3g (%.3g of truth) rhat %.2f sdgp rhat %.2f\n", s, median(b[, , 1]), median(b[, , 1]) / s, posterior::rhat(b[, , 1]), posterior::rhat(b[, , 2])))
  }
  flush.console()
}
cat("DONE\n")
