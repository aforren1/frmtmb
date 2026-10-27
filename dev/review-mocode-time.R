## Fix round 1 check: does mo_check_fit() (logLik(), summary()) or
## lp_delta_mo() (fitted() standard errors) cost time on a large new
## fit? 0-based codes so the base arm fits the same model. Seed 9,
## n = 1e5. Minimum of 5 rounds per block, with a fixed arithmetic
## control that must read the same in both arms.
##   Rscript dev/review-mocode-time.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", arm, "lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(9)
n <- 1e5
d <- data.frame(k0 = sample(0:3, n, TRUE), j0 = sample(0:4, n, TRUE),
                z = rnorm(n))
d$y <- 1 + c(0, 1, 1.5, 1.7)[d$k0 + 1] + c(0, 0.3, 0.5, 0.6, 0.65)[d$j0 + 1] +
  0.5 * d$z + rnorm(n)
fit <- suppressWarnings(suppressMessages(frm(y ~ mo(k0) + mo(j0) + z,
                                             data = d)))
tm <- function(expr, reps) {
  f <- function() { for (i in seq_len(reps)) force(expr()) }
  min(vapply(1:5, function(r) system.time(f())[["elapsed"]], 0))
}
ctrl <- tm(function() sum(sort(runif(2e5))), 200)
ll <- tm(function() logLik(fit), 50000)
sm <- tm(function() summary(fit), 200)
fx <- tm(function() fixef(fit), 2000)
ft <- tm(function() fitted(fit), 60)
cat(sprintf("control %.3f s | logLik x50000 %.3f s | summary x200 %.3f s | fixef x2000 %.3f s | fitted (SE, n = 1e5) x60 %.3f s\n",
            ctrl, ll, sm, fx, ft))
cat("fitted SE head:", format(head(fitted(fit)[, 2], 3), digits = 7), "\n")
