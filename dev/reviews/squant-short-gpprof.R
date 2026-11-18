# Short check of lane squant fix round 1, item 4 follow-up: where the
# 129 s of the quantile_elf() gp(x) + x2 fit on 40 rows go
# (squant-short-gp.txt; smoothing 1 takes 0.6 s). Rprof, top self and
# total times.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-gpprof.R > dev/reviews/squant-short-gpprof.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261117)
n <- 400
d <- data.frame(x = runif(n), x2 = runif(n),
                g = factor(sample(sprintf("g%02d", 1:25), n, TRUE)))
re0 <- rnorm(25, 0, 0.6)
re1 <- rnorm(25, 0, 0.5)
gi <- as.integer(d$g)
d$y <- 1 + 2 * d$x + re0[gi] + re1[gi] * d$x + (0.5 + d$x) * rnorm(n)
d$ys <- sin(2 * pi * d$x) + 0.5 * d$x2^2 + 0.4 * rnorm(n)
pf <- "dev/reviews/squant-short-gpprof.out"
Rprof(pf, interval = 0.05)
t0 <- proc.time()[["elapsed"]]
fit <- withCallingHandlers(
  frm(bf(ys ~ gp(x) + x2, quantile = 0.5), family = quantile_elf(),
      data = d[1:40, ]),
  warning = function(w) {
    cat("W:", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  })
Rprof(NULL)
cat(sprintf("elapsed %.1f s, conv %d, logLik %.4f, fn evals %s, iterations %s\n",
            proc.time()[["elapsed"]] - t0, fit$opt$convergence,
            as.numeric(logLik(fit)),
            paste(fit$opt$evaluations, collapse = "/"),
            if (is.null(fit$opt$iterations)) "NA" else fit$opt$iterations))
s <- summaryRprof(pf)
cat("== by.total, top 30 ==\n")
print(head(s$by.total, 30))
cat("== by.self, top 15 ==\n")
print(head(s$by.self, 15))
unlink(pf)
cat("DONE\n")
