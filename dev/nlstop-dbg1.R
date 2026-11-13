# Lane nlstop: debug of the first build on items 2584 and 4705.
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
# 2584
set.seed(8)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, verbose = FALSE))
d$gr <- factor(rep(1:20, length.out = nrow(d)))
set.seed(8)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, verbose = FALSE))
d$z <- runif(nrow(d)); d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
d$g8 <- factor(sample(letters[1:8], nrow(d), TRUE))
d$gr <- factor(rep(1:20, length.out = nrow(d)))
d$y <- d$y + 0.6 * rnorm(20)[d$gr]
f <- frm(bf(y ~ s(x1), sigma ~ s(x2)), data = d)
cat("random:", length(f$obj$env$random), " sd index:",
    tryCatch(ns$log_sd_theta_index(f), error = conditionMessage), "\n")
om <- ns$outer_par_map(f)
print(om)
for (k in seq_along(f$estimates$theta)) {
  cat("k", k, "scale", ns$se_sd_scale(f, k), "\n")
}
p <- f$opt$par
for (j in which(om$comp == "theta")) {
  q <- p; q[j] <- q[j] - 10
  cat("theta pos", j, "f(p)", f$opt$objective, "f(q)", f$obj$fn(q), "\n")
}
