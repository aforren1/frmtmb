# Lane illse fix round 1, item 7: time per function of a wide fit (the
# reviewer's design at a smaller n), inclusive, by trace().
#   Rscript dev/illse-fx-cost2.R base|lane|dev <re: 0|1> <n>
args <- commandArgs(TRUE)
arm <- args[1]
re <- as.integer(args[2])
n <- as.numeric(args[3])
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-illse-lib"),
            file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
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
fo <- stats::as.formula(paste("y ~", paste(c(paste0("x", 1:6), colnames(Z)),
                                          collapse = " + "),
                              if (re) "+ (1 | g)"))
tm <- new.env()
fns <- c("newton_polish", "conv_headroom", "grad_headroom", "grad_verdict",
         "check_convergence", "se_check", "cov_from_hessian",
         "se_reparam_cov", "sdr_of", "fit_end_checks", "escape_stationary",
         "fd_hessian_noise", "se_tier3")
fns <- fns[vapply(fns, exists, TRUE, envir = ns)]
for (fn in fns) {
  local({
    f0 <- fn
    trace(f0, where = ns, print = FALSE,
          tracer = bquote(assign(paste0(".t0_", .(f0)),
                                 c(get0(paste0(".t0_", .(f0)), envir = tm, inherits = FALSE,
                                        ifnotfound = NULL),
                                   proc.time()[["elapsed"]]), envir = tm)),
          exit = bquote({
            st <- get(paste0(".t0_", .(f0)), envir = tm)
            el <- proc.time()[["elapsed"]] - st[length(st)]
            assign(paste0(".t0_", .(f0)), st[-length(st)], envir = tm)
            if (!length(st[-length(st)])) {
              assign(.(f0), c(get0(.(f0), envir = tm, inherits = FALSE, ifnotfound = NULL),
                              el), envir = tm)
            }
          }))
  })
}
t0 <- proc.time()[["elapsed"]]
f <- suppressWarnings(frm(fo, data = d))
cat(sprintf("arm %s | re %d | n %g | frm %.1f s | code %d\n", arm, re, n,
            proc.time()[["elapsed"]] - t0, f$opt$convergence))
for (fn in fns) {
  x <- get0(fn, envir = tm, inherits = FALSE, ifnotfound = NULL)
  if (length(x)) cat(sprintf("  %-20s calls %d total %.1f s\n", fn,
                             length(x), sum(x)))
}
