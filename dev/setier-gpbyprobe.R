# Lane setier: test-gp-by.R's non-isotropic gp(x, z, by = f) fit: which
# hyperparameters are flat, their values and blocks, and how the
# objective moves when each is pushed either way.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                w = stats::runif(n, 0.5, 2))
d$y <- 0.5 + ifelse(d$f == "a", sin(d$x),
                    ifelse(d$f == "b", cos(d$x), 0.2 * d$x)) +
  stats::rnorm(n, 0, 0.3)
d$z <- stats::runif(nrow(d), 0, 3)
f2 <- suppressWarnings(
  frm(bf(y ~ gp(x, z, by = f, iso = FALSE, k = 5)), data = d))
print(round(f2$estimates$theta, 3))
for (bk in f2$frame$re_blocks) {
  cat(bk$term_label, bk$covstruct, "theta", bk$theta_idx, "\n")
}
print(ns$sdr_of(f2)$se_lost)
p <- f2$opt$par
f0 <- f2$obj$fn(p)
nm <- ns$outer_par_names(f2)
for (k in c("theta_7", "theta_9")) {
  j <- match(k, nm)
  for (s in c(-2, 2)) {
    q <- p
    q[j] <- q[j] + s
    cat(k, s, f2$obj$fn(q) - f0, "\n")
  }
}
