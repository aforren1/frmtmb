# Lane gpfield fix round: Rprof of predict(no newdata, ndraws = 100) on
# dev/gpfield-cost.R's poisson n 300 fit. Args: arm
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(a[1])
n <- 300
set.seed(9300 + n)
d <- data.frame(x = sort(round(stats::runif(n, 0, 50), 4)))
d$y <- stats::rpois(n, exp(1 + sin(d$x / 3)))
fit <- gf_fit(bf(y ~ gp(x)), data = d, family = poisson())$fit
set.seed(1)
invisible(predict(fit, summary = FALSE, ndraws = 10))
pf <- file.path(GF_OUT, paste0("predict-prof-", a[1], ".out"))
Rprof(pf, interval = 0.005)
t0 <- proc.time()[[3]]
set.seed(1)
p <- predict(fit, summary = FALSE, ndraws = 100)
cat("PREDPROF", a[1], "100 draws", proc.time()[[3]] - t0, "s\n")
Rprof(NULL)
s <- summaryRprof(pf)
print(head(s$by.total, 30))
print(head(s$by.self, 15))
