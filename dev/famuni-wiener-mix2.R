# Where the logLik of 7713.9 in dev/famuni-out/wiener-mix.txt comes from:
# the mixture(wiener(link_ndt = "log", allow_unreachable = TRUE),
# lognormal()) fit on data with no fast guesses, seed 20261013.
# Output: dev/famuni-out/wiener-mix2.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
seed <- 20261013L
set.seed(seed)
n <- 400
d <- ddm_simulate(n, mu = 0.8, bs = 1.3, ndt = 0.3, bias = 0.45)
f <- suppressWarnings(frm(bf(rt | dec(upper) ~ 1),
         family = frmtmb::mixture(wiener(link_ndt = "log",
                                         allow_unreachable = TRUE),
                                  lognormal()), data = d))
print(f$opt$par)
dp <- c("mu1", "bs1", "ndt1", "bias1", "mu2", "sigma2", "theta1", "theta2")
v <- lapply(dp, function(p) {
  r <- tryCatch(frm_linpred(f, dpar = p, type = "response")[1],
                error = function(e) NA_real_)
  r
})
names(v) <- dp
str(v)
w1 <- wiener(link_ndt = "log", allow_unreachable = TRUE)
l1 <- w1$lpdf(d$rt, list(mu = v$mu1, bs = v$bs1, ndt = v$ndt1,
                         bias = v$bias1), list(dec = d$upper))
cat("wiener component lpdf at the optimum: range", range(l1), "\n")
l2 <- dlnorm(d$rt, v$mu2, v$sigma2, log = TRUE)
cat("lognormal component lpdf: range", range(l2), " sum", sum(l2), "\n")
cat("objective at the optimum", f$obj$fn(f$opt$par), "\n")
