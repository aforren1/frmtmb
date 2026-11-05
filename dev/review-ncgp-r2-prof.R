# Fix round 1, item 8: where the 33 s go at 500 positions, profiled. Rprof by.total and counts of the
# tape's calls, one fit, in the fit's own stages. Seed 20261103, as dev/ncgp-perf.R. Args: n
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
n <- as.integer(args[1])
ns <- asNamespace("frmtmb")
set.seed(20261103); npos <- min(n, 500L); pos <- sort(runif(npos, 0, 10)); x <- rep(pos, length.out = n)
d <- data.frame(x = x, y = sin(x) + 0.3 * rnorm(n))
pf <- file.path(REV_OUT, sprintf("r2-prof-%d.out", n))
t0 <- proc.time()[[3]]
Rprof(pf, interval = 0.05, memory.profiling = FALSE)
fit <- suppressWarnings(frm(y ~ gp(x), d, control = frmtmb_control(verbose = TRUE)))
Rprof(NULL)
cat(sprintf("RPF total %.1f s ll %.4f iter %d evals %s\n", proc.time()[[3]] - t0, logLik(fit), fit$opt$iterations, paste(fit$opt$evaluations, collapse = "/")))
s <- summaryRprof(pf)$by.total
keep <- c("gp_marginal_max", "gp_marginal_opt", "gp_marginal_polish", "gp_start_scan", "gp_tape_scope", "optimize_obj", "nlminb", "RTMB::MakeADFun", "MakeADFun", "autoscale_sdreport", "RTMB::sdreport", "sdreport",
          "check_convergence", "conv_headroom", "fit_outer_hessian", "fit_outer_hessian_at", "se_check", "se_analysis", "se_report", "se_tier3", "gp_nugget_message", "gp_position_info", "newton_polish", "optimHess", "obj$fn", "obj$gr", "fn", "gr", "f", "EvalADFunObject", "chol", "post_fit", "fit_assembled", "assemble_frame", "frm")
print(round(s[rownames(s) %in% paste0('"', keep, '"'), c("total.time", "total.pct")], 1))
print(utils::head(round(s[, c("total.time", "total.pct")], 1), 60))
# the tape's costs, on the finished fit
p <- fit$opt$par
tm <- function(e) { t0 <- proc.time()[[3]]; force(e); proc.time()[[3]] - t0 }
cat(sprintf("RPF one fn %.2f gr %.2f sdreport %.2f\n", tm(fit$obj$fn(p + 1e-5)), tm(fit$obj$gr(p + 2e-5)), tm(RTMB::sdreport(fit$obj))))
# the marginal's exact Hessian against the outer differences, for the cut
fr <- fit$frame
s0 <- ns$gp_marginal_setup(fr)
est <- fit$estimates[c("beta", "betad", "theta")]
old <- RTMB::TapeConfig(); RTMB::TapeConfig(matmul = "atomic")
t1 <- tm(ob <- RTMB::MakeADFun(ns$gp_marginal_nll(s0), list(beta = est$beta, betad = est$betad, theta = est$theta), silent = TRUE))
t2 <- tm(H <- optimHess(ob$par, ob$fn, ob$gr))
RTMB::TapeConfig(old)
cat(sprintf("RPF marginal tape %.2f s, Hessian by its gradients %.2f s, value at the optimum %.8f vs fit %.8f\n", t1, t2, ob$fn(ob$par), fit$opt$objective))
Vm <- solve(H)
Vz <- as.matrix(vcov(fit, full = TRUE))
cat("RPF SE from the marginal", signif(sqrt(diag(Vm)), 5), "\nRPF SE from the fit     ", signif(sqrt(diag(Vz)), 5), "\n")
cat("DONE\n")
