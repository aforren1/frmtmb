# Review of lane adhess, item 4: cost of the AD Hessian at p = 200 and
# p = 500 without random effects, against 0.71.0's differences.
#   REV_ARM=lane|base Rscript dev/adhess-rev-cost.R <p> <n> [ignore]
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
a <- commandArgs(TRUE)
p <- as.integer(a[1])
n <- as.integer(a[2])
ign <- isTRUE(a[3] == "ignore")
peak <- function() ps::ps_memory_info()[["peak_wset"]] / 2^20
set.seed(100 + p)
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- sprintf("x%03d", seq_len(p))
b <- rnorm(p, 0, 0.3 / sqrt(p))
d <- as.data.frame(X)
d$y <- rbinom(n, 1, plogis(0.2 + X %*% b))
f <- stats::as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
cat(sprintf("p %d n %d check_se %s | peak before fit %.0f MB\n", p, n,
            if (ign) "ignore" else "default", peak()))
ctl <- if (ign) frmtmb_control(check_se = "ignore") else frmtmb_control()
t_fit <- system.time(fit <- suppressWarnings(
  frm(f, family = bernoulli(), data = d, control = ctl)))[["elapsed"]]
same <- identical(unname(fit$opt$par), unname(ns$sdr_outer_point(fit$obj)))
cat(sprintf("fit %.2f s | code %d | opt$par == last.par.best %s | peak %.0f MB\n",
            t_fit, fit$opt$convergence, same, peak()))
t_v <- system.time(V <- vcov(fit))[["elapsed"]]
cat(sprintf("first vcov %.2f s | peak %.0f MB\n", t_v, peak()))
at <- ns$sdr_outer_point(fit$obj)
st <- ns$obj_state_save(fit$obj)
t_he <- system.time(Hx <- fit$obj$he(at))[["elapsed"]]
cat(sprintf("obj$he alone %.2f s | peak %.0f MB\n", t_he, peak()))
t_fd <- system.time(Hf <- ns$fd_hessian(at, fit$obj$gr)$H)[["elapsed"]]
ns$obj_state_restore(fit$obj, st)
cat(sprintf("fd_hessian alone %.2f s (2p = %d gradients) | peak %.0f MB\n",
            t_fd, 2L * p, peak()))
D <- sqrt(diag(Hx))
cat(sprintf("max |H_ad - H_fd| / D D' %.2e\n",
            max(abs((Hx + t(Hx)) / 2 - Hf) / outer(D, D))))
