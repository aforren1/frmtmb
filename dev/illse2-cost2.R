# Lane illse2, item 5: every Hessian the wide design of dev/illse2-cost.R
# takes, with its time, its caller and the point (a hash), so that the
# repeated ones show. reparam_hessian() is one pass in the design's
# coordinates (4 k gradients with Richardson); fd_hessian() is the
# natural finite-difference one; obj$he() under obj_he_cached().
#   ILLSE2_ARM=lane|base|r8 Rscript dev/illse2-cost2.R <re 0|1> <n>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
re <- as.integer(a[1])
n <- as.numeric(a[2])
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
log <- list()
callers <- c("newton_polish", "conv_headroom", "grad_headroom",
             "se_reparam_cov", "check_convergence", "sd_up_escape",
             "headroom_reparam_hessian", "se_check", "cov_from_hessian",
             "fit_outer_hessian", "sdr_rescue_hessian", "se_tier3")
for (fn in intersect(c("reparam_hessian", "fd_hessian", "optimHess"),
                     ls(ns, all.names = TRUE))) {
  local({
    f0 <- fn
    trace(f0, where = ns, print = FALSE,
          tracer = quote(.i2_t0 <- proc.time()[["elapsed"]]),
          exit = bquote({
            who <- vapply(sys.calls(), function(cl) {
              h <- cl[[1L]]
              if (is.name(h)) as.character(h) else ""
            }, "")
            pp <- if (exists("p", inherits = FALSE)) p else NA
            assign("log", c(get("log", envir = globalenv()), list(list(
              fn = .(f0), t = proc.time()[["elapsed"]] - .i2_t0,
              who = paste(who[who %in% callers], collapse = ">"),
              p = sum(as.numeric(pp) * seq_along(pp))))),
              envir = globalenv())
          }))
  })
}
t0 <- proc.time()[["elapsed"]]
r <- cap(frm(fo, data = d))
t1 <- proc.time()[["elapsed"]] - t0
f <- r$v
cat(sprintf(paste0("arm %s | re %d | n %g | outer %d | frm %.1f s | ",
                   "code %d | lost %d | root %s | warn %s | evals %s\n"),
            arm, re, n, length(f$opt$par), t1, f$opt$convergence,
            length(ns$sdr_of(f)$se_lost),
            !is.null(ns$sdr_of(f)$cov_fixed_root), wtag(r$w),
            format(f$opt$evals)))
for (e in log) {
  cat(sprintf("  %-16s %6.1f s | point %.12g | %s\n", e$fn, e$t, e$p, e$who))
}
cat(sprintf("  total in traced Hessians %.1f s\n",
            sum(vapply(log, `[[`, 0, "t"))))
