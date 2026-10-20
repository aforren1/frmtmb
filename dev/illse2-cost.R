# Lane illse2, item 5: where a wide ill-conditioned fit spends its time.
# The reviewer's design of dev/illse-rev-cost.R: n rows, 100 outer
# parameters (intercept, a raw polynomial of degree 6 on [1, 2], 92
# standard normal covariates, sigma), or with (1 | g) over 100 groups
# (91 covariates and the sd). One frm() under Rprof (interval 0.02 s);
# prints the wall time, the inclusive time of the functions the SE and
# convergence machinery runs, and how many gradients each phase asked.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-cost.R <re 0|1> <n> [prof 0|1]
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
# the difference scheme of the experiment (item 5); unset, the default
ri <- Sys.getenv("I2_RICH", "")
if (nzchar(ri)) options(frmtmb.i2_rich = if (ri == "single") "single" else ri == "1")
a <- commandArgs(TRUE)
re <- as.integer(a[1])
n <- as.numeric(a[2])
prof <- if (length(a) > 2) as.integer(a[3]) == 1L else TRUE
set.seed(1)
nc <- if (re) 91L else 92L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(sample(1:100, n, TRUE)))
for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
Z <- matrix(stats::rnorm(n * nc), n, nc)
colnames(Z) <- paste0("z", seq_len(nc))
d <- cbind(d, Z)
d$y <- sin(3 * d$x) + 0.05 * rowSums(Z[, 1:10]) +
  (if (re) stats::rnorm(100, 0, 0.3)[d$g] else 0) + stats::rnorm(n, 0, 0.5)
fo <- stats::as.formula(paste("y ~", paste(c(paste0("x", 1:6),
                                             colnames(Z)), collapse = " + "),
                              if (re) "+ (1 | g)"))
pf <- file.path("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-log",
                sprintf("cost-%s-%d-%g.rprof", arm, re, n))
if (prof) Rprof(pf, interval = 0.02)
pt0 <- proc.time()
t0 <- pt0[["elapsed"]]
r <- cap(frm(fo, data = d))
pt1 <- proc.time() - pt0
t1 <- pt1[["elapsed"]]
cpu <- pt1[["user.self"]] + pt1[["sys.self"]]
if (prof) Rprof(NULL)
f <- r$v
cat(sprintf(paste0("arm %s | re %d | n %g | outer %d | frm %.1f s (cpu %.1f) | ",
                   "code %d | lost %d | root %s | warn %s | evals %s\n"),
            arm, re, n, length(f$opt$par), t1, cpu, f$opt$convergence,
            length(ns$sdr_of(f)$se_lost),
            !is.null(ns$sdr_of(f)$cov_fixed_root), wtag(r$w),
            format(f$opt$evals)))
if (prof) {
  s <- summaryRprof(pf)$by.total
  keep <- c("frm", "fit_assembled", "optimize_obj", "run_optimizer",
            "newton_polish", "sd_up_escape", "check_convergence",
            "conv_headroom", "grad_headroom", "headroom_drop_rows",
            "se_sd_gain_up", "fit_outer_hessian", "fit_outer_hessian_at",
            "obj_he_cached", "reparam_hessian_refined_at",
            "reparam_hessian", "se_check", "se_reparam_cov",
            "cov_from_hessian", "sdr_rescue", "autoscale_sdreport",
            "sdreport", "fd_hessian", "se_tier3", "se_obj_curv",
            "assemble_frame", "MakeADFun", "predictor_resid_sd",
            "start_ls_plan", "optimHess", "eigen", "qr", "chol",
            "nlminb")
  hit <- s[rownames(s) %in% paste0("\"", keep, "\""), , drop = FALSE]
  hit <- hit[order(-hit$total.time), c("total.time", "total.pct")]
  print(hit)
}
