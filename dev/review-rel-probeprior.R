# Release review 0.80.0: merge step 3 on the multistart message. Under a
# prior with a density the message must say "penalized log-likelihood"
# and print the maximized values; without one, "log-likelihood". The
# NEWS case `b1 * (1 - exp(-b2 * x))` under gaussian(), seeds 1 to 5.
#
#   Rscript dev/review-rel-probeprior.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
fm <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
pr <- set_prior("normal(0, 5)", nlpar = "b1") +
  set_prior("normal(0, 5)", nlpar = "b2")
for (seed in 1:5) {
  set.seed(seed)
  d <- data.frame(x = runif(200, 0, 4))
  d$y <- 3 * (1 - exp(-1.2 * d$x)) + rnorm(200, 0, 0.3)
  for (arm in c("no prior", "prior")) {
    ms <- character()
    f <- withCallingHandlers(
      if (arm == "prior") frm(fm, data = d, prior = pr) else frm(fm, data = d),
      message = function(m) {
        ms <<- c(ms, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        ms <<- c(ms, paste("WARNING", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    hit <- grep("sign|restart", ms, value = TRUE)
    got <- regmatches(hit, gregexpr("log-likelihood -?[0-9.]+", hit))
    nums <- as.numeric(sub(".*log-likelihood ", "", unlist(got)))
    cat(sprintf("seed %d %-8s pen %.4f logLik %.4f | words: %s | last printed %s\n",
                seed, arm, as.numeric(logLik_penalized(f)),
                as.numeric(logLik(f)),
                paste(unique(unlist(regmatches(hit, gregexpr(
                  "(penalized )?log-likelihood", hit)))), collapse = ","),
                if (length(nums)) format(nums[length(nums)]) else "-"))
  }
}
cat("PROBEPRIOR DONE\n")
