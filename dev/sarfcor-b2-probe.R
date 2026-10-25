# Part B2: the weights check and a registered full-rank structure. The
# model of test-weights-unbounded.R's "a registered structure is checked
# by its own rank": Gamma(link = "log"), weights 1.5 on every row, and an
# observation-level block, built in, (1 | obs), and registered, myiid(1
# | obs), the same iid block under another name. wu_data(), seed
# 20261020, 72 rows. Run once per library:
#   Rscript dev/sarfcor-b2-probe.R base > dev/sarfcor-out/b2-base.txt
#   Rscript dev/sarfcor-b2-probe.R lane > dev/sarfcor-out/b2-lane.txt
arm <- commandArgs(TRUE)[1]
lib <- c(if (identical(arm, "lane")) "C:/Users/adf44/source/r/wt-sarfcor-lib",
         "C:/Users/adf44/source/r/rellib-r12",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressMessages(library(frmtmb))
cat("arm:", arm, "lib:", find.package("frmtmb"), "\n")
frmtmb_register_covstruct(
  "myiid",
  lpdf = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
  },
  vcov = function(theta, blk) {
    matrix(exp(2 * theta[1]), 1, 1,
           dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  },
  npar = 1L, sd_idx = 1L)
set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                tt = factor(rep(seq_len(nt), ng)),
                x = rnorm(n), xs = rnorm(n), w = runif(n, 0.5, 2))
d$obs <- factor(seq_len(n))
eta <- 1 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$w15 <- 1.5
G <- Gamma(link = "log")
for (f in list(y | weights(w15) ~ x + (1 | obs),
               y | weights(w15) ~ x + myiid(1 | obs))) {
  w <- character()
  r <- tryCatch(withCallingHandlers(frm(f, family = G, data = d),
                                    warning = function(x) {
                                      w <<- c(w, conditionMessage(x))
                                      invokeRestart("muffleWarning")
                                    }),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(deparse(f), "-> REFUSED:", substr(conditionMessage(r), 1, 110), "\n")
  } else {
    cat(deparse(f), "-> FIT logLik", format(as.numeric(logLik(r)),
                                           digits = 10),
        "shape", format(exp(r$estimates$betad[1]), digits = 4),
        "warnings:", length(w), "\n")
    if (length(w)) cat("  first warning:", substr(w[1], 1, 120), "\n")
  }
}
