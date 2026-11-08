# Short check of lane latperf's fix round, m3: the saved fit after the readers
# holds no selected-inverse lookup, and the readers of a fit read back with
# readRDS() give identical values. Design "two" of dev/latperf-prof.R, seed 1.
# Usage: Rscript dev/review-latperf-saved.R <n>
n <- as.integer(commandArgs(TRUE)[1])
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat(sprintf("lib: %s frmtmb %s\n", dirname(find.package("frmtmb")),
            format(packageVersion("frmtmb"))))
mb <- function(x) length(serialize(x, NULL)) / 2^20
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n),
                id = factor(sample.int(50, n, TRUE)))
e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
u <- rnorm(50, 0, 0.3)
d$y <- rpois(n, exp(.5 + .3 * d$x + e + u[d$id]))
fit <- suppressWarnings(frm(y ~ x + (1 | id) + ar(t, g, cov = TRUE),
                            family = poisson(), data = d))
before <- mb(fit)
set.seed(2)
nd <- data.frame(t = n + seq_len(50), g = factor(1), x = rnorm(50),
                 id = factor(sample.int(50, 50, TRUE), levels = 1:50))
rd <- function(f) {
  list(fitted = suppressWarnings(fitted(f)),
       resid = suppressWarnings(residuals(f)),
       linpred = frm_linpred(f, se.fit = TRUE),
       newdata = suppressWarnings(fitted(f, newdata = nd)))
}
t1 <- system.time(a <- rd(fit))[["elapsed"]]
S <- fit$cache$Vjoint$V
cat(sprintf("RESULT n %d readers %.2f s; saved before %.2f after %.2f MB\n",
            n, t1, before, mb(fit)))
cat(sprintf("RESULT Vjoint %.2f MB: elements %s\n", mb(fit$cache$Vjoint),
            paste(names(unclass(S)), collapse = ",")))
cat("RESULT any environment in V:", any(vapply(unclass(S), is.environment,
                                               TRUE)), "\n")
tf <- tempfile(fileext = ".rds")
saveRDS(fit, tf)
cat(sprintf("RESULT saveRDS file %.2f MB\n", file.size(tf) / 2^20))
f2 <- readRDS(tf)
cat("RESULT reloaded sparse:", frmtmb:::jc_is_sparse(f2$cache$Vjoint), "\n")
t2 <- system.time(b <- rd(f2))[["elapsed"]]
for (k in names(a)) {
  cat(sprintf("RESULT %-8s identical after readRDS: %s\n", k,
              identical(a[[k]], b[[k]])))
}
cat(sprintf("RESULT reloaded readers %.2f s\n", t2))
unlink(tf)
