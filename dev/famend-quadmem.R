# Lane famend, 3580: the time and peak working set of a censored poisson
# (1 | g) fit under Laplace and single-block quadrature, the release
# review's construction (dev/rel073-rev-quadmem.R, seed 20261021, 6 rows
# per group, right-censored at 5), on one build. One fit per process, so
# the peak is that fit's.
#   Rscript dev/famend-quadmem.R <lib or "base"> <n> <quadrature TRUE|FALSE>
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressMessages(library(frmtmb))
n <- as.integer(a[2])
quad <- as.logical(a[3])
set.seed(20261021)
ng <- n / 6
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = 6)))
u <- rnorm(ng, 0, 0.6)[d$g]
d$y <- rpois(n, exp(1 + 0.3 * d$x + u))
d$cc <- ifelse(d$y >= 5, "right", "none")
d$y <- pmin(d$y, 5)
gc()
m0 <- ps::ps_memory_info(ps::ps_handle())[["peak_wset"]]
t0 <- proc.time()[[3]]
f <- suppressWarnings(frm(bf(y | cens(cc) ~ x + (1 | g)), family = poisson(),
                          data = d, quadrature = quad))
t1 <- proc.time()[[3]]
m1 <- ps::ps_memory_info(ps::ps_handle())[["peak_wset"]]
cat(sprintf(paste0("QUADMEM frmtmb %s n %d quadrature %s: logLik %.6f | ",
                   "%.1f s | peak working set %.0f MB (%.0f MB before)\n"),
            format(packageVersion("frmtmb")), n, quad,
            as.numeric(logLik(f)), t1 - t0, m1 / 2^20, m0 / 2^20))
