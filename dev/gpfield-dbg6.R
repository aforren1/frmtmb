# Lane gpfield: test-brms-lik-gp-exact.R's approximate row, y ~ gp(x, k =
# 10) on gpx_data() (set.seed(20261102)): the default fit and the profile
# of the log-likelihood over the length scale (each point the best of the
# other parameters from two sd starts). Args: arm kernel
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(a[1])
cv <- if (length(a) > 1) a[2] else "exp_quad"
set.seed(20261102)
n <- 80
d <- data.frame(x = sample(0:12, n, TRUE), z = sample(0:12, n, TRUE),
                f = factor(sample(c("a", "b"), n, TRUE)))
e <- stats::rnorm(n, 0, 0.3)
d$y <- 2 * sin(0.7 * d$x) + e
fml <- bf(as.formula(bquote(y ~ gp(x, k = 10, cov = .(cv)))))
q <- gf_fit(fml, data = d, family = gaussian())
f <- q$fit
bk <- f$frame$re_blocks[[1]]
cat("DEFAULT ll", as.numeric(logLik(f)), "theta", f$estimates$theta,
    "spacing", 2 * 1.25 / 10,
    "| warn", sub(":.*", "", q$warn), "\n")
th <- which(names(f$opt$par) == "theta")
for (ls in c(1e-4, 1e-3, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 0.8)) {
  best <- -Inf
  for (s0 in c(0.5, 5)) {
    p0 <- f$opt$par
    p0[th] <- c(log(s0), log(ls))
    fn <- function(q) { p <- p0; p[-th[2]] <- q; f$obj$fn(p) }
    gr <- function(q) { p <- p0; p[-th[2]] <- q; f$obj$gr(p)[-th[2]] }
    o <- tryCatch(nlminb(p0[-th[2]], fn, gr), error = function(e) NULL)
    if (!is.null(o)) best <- max(best, -o$objective)
  }
  cat(sprintf("PROF ls %.4g (x %.3g of the spacing) ll %.4f\n", ls,
              ls / 0.25, best))
}
