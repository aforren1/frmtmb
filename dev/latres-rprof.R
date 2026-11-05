# Where a frm() call on a latent residual spends its time. Rprof over one fit,
# with the top self and total entries. Seed 167.
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
args <- commandArgs(TRUE)
G <- if (length(args)) as.integer(args[[1]]) else 200L
k <- if (length(args) > 1L) as.integer(args[[2]]) else 5L
set.seed(167)
d <- expand.grid(t = seq_len(k), g = factor(seq_len(G)))
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), k, sd = 0.5)))
d$y <- rpois(n, exp(0.5 + 0.3 * d$x + e))
f <- tempfile()
Rprof(f, interval = 0.05)
t0 <- proc.time()[[3]]
fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(), data = d)
print(proc.time()[[3]] - t0)
Rprof(NULL)
s <- summaryRprof(f)
print(utils::head(s$by.self, 15))
tt <- s$by.total
print(utils::head(tt[grepl("frm|se_|sdr|check|newton|hess|outer|fit_|frame", rownames(tt)), ], 40))
