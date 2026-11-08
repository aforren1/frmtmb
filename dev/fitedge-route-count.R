# Item 3: how often a fit asks for the marginal's gradient, and what one ask
# costs, at `positions` (data as dev/fitedge-hess-time.R). Args: arm, positions.
a <- commandArgs(TRUE)
arm <- a[1]; npos <- as.integer(a[2])
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
orig <- get("gp_marginal_outer_gr", envir = ns)
calls <- 0L
tsum <- 0
unlockBinding("gp_marginal_outer_gr", ns)
assign("gp_marginal_outer_gr", function(fit, at) {
  calls <<- calls + 1L
  t0 <- proc.time()[[3]]
  on.exit(tsum <<- tsum + proc.time()[[3]] - t0)
  orig(fit, at)
}, envir = ns)
lockBinding("gp_marginal_outer_gr", ns)
set.seed(20261103)
pos <- sort(runif(npos, 0, 10))
d <- data.frame(x = pos, y = sin(pos) + 0.3 * rnorm(npos))
t0 <- proc.time()[[3]]
f <- suppressWarnings(frm(y ~ gp(x), data = d))
t1 <- proc.time()[[3]]
cat(sprintf("RC pos %d frm %.2f s: %d asks, %.2f s inside them\n", npos, t1 - t0, calls, tsum))
s <- summary(f)
cat(sprintf("RC after summary(): %d asks, %.2f s\n", calls, tsum))
Rprof(file.path(FE_OUT, sprintf("route-prof-%d.out", npos)), interval = 0.01)
g <- orig(f, f$opt$par)
Rprof(NULL)
p <- summaryRprof(file.path(FE_OUT, sprintf("route-prof-%d.out", npos)))$by.total
print(head(round(p[, c("total.time", "total.pct")], 2), 25))
cat("DONE\n")
