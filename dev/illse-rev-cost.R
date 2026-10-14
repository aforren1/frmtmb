# Reviewer of lane illse, attack 8: the cost of the design-coordinate
# Hessian on a wide design. n = 1e5, 100 outer parameters (intercept, a
# raw polynomial of degree 6 on [1, 2], 92 standard normal covariates,
# sigma), with and without (1 | g) (99 parameters plus the sd then: 91
# covariates). Times frm() and, on the lane, the parts the lane added.
#   REV_ARM=lane|base Rscript dev/illse-rev-cost.R <re: 0|1>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
re <- as.integer(commandArgs(TRUE)[1])
set.seed(1)
n <- 1e5
nc <- if (re) 91L else 92L
d <- data.frame(x = stats::runif(n, 1, 2), g = factor(sample(1:100, n, TRUE)))
for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
Z <- matrix(stats::rnorm(n * nc), n, nc)
colnames(Z) <- paste0("z", seq_len(nc))
d <- cbind(d, Z)
d$y <- sin(3 * d$x) + 0.05 * rowSums(Z[, 1:10]) +
  (if (re) stats::rnorm(100, 0, 0.3)[d$g] else 0) + stats::rnorm(n, 0, 0.5)
fo <- stats::as.formula(paste("y ~", paste(c(paste0("x", 1:6), colnames(Z)),
                                          collapse = " + "),
                              if (re) "+ (1 | g)"))
tm <- list()
if (rev_arm == "lane") {
  for (fn in c("se_reparam_cov", "reparam_hessian_refined", "conv_headroom",
               "newton_polish")) {
    local({
      f0 <- fn
      trace(f0, tracer = bquote(assign(".t0", proc.time()[["elapsed"]],
                                       envir = globalenv())),
            exit = bquote(tm[[.(f0)]] <<- c(tm[[.(f0)]],
                                           proc.time()[["elapsed"]] -
                                             get(".t0", envir = globalenv()))),
            print = FALSE, where = ns)
    })
  }
}
t0 <- proc.time()[["elapsed"]]
r <- cap(frm(fo, data = d))
t1 <- proc.time()[["elapsed"]] - t0
f <- r$v
root <- ns$sdr_of(f)$cov_fixed_root
cat(sprintf("re %d | n %d | outer %d | frm %.1f s | code %d | lost %d | root %s | warn %s\n",
            re, n, length(f$opt$par), t1, f$opt$convergence,
            length(ns$sdr_of(f)$se_lost), !is.null(root), wtag(r$w)))
for (nm in names(tm)) cat(sprintf("  %-24s calls %d total %.1f s\n", nm,
                                  length(tm[[nm]]), sum(tm[[nm]])))
t2 <- proc.time()[["elapsed"]]
v <- suppressWarnings(fitted(f, newdata = d[1:1000, ], re_formula = NA))
cat(sprintf("  fitted(1000 rows) %.2f s\n", proc.time()[["elapsed"]] - t2))
