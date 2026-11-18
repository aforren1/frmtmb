# Lane sparsese: the sparse route against the dense one on one fit.
# A gaussian REML fit, degree-7 raw polynomial, (1 | g) with `ng` groups of
# 2 rows, fitted twice in one process: once with
# options(frmtmb.joint_root_max = 1e6) (the dense factor, sdr_unwhiten() and
# joint_cov_white_inner()) and once with 100 (sparse solves in the whitened
# coordinates). Prints, for each reader, the largest relative difference
# between the two, and the time of each reader's first call.
#   Rscript sparsese-agree.R <groups> <seed> [group sd, default 0]
a <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(a[1])
seed <- as.integer(a[2])
gsd <- if (length(a) >= 3L) as.numeric(a[3]) else 0
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "groups", ng, "seed", seed, "\n")
quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (!grepl("a product such as a' V a", m, fixed = TRUE)) {
      cat("   W:", substr(m, 1, 160), "\n")
    }
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
}
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, gsd)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
# the sdreport is memoized on the fit at its first use, so each arm's limit
# is in force while it is formed and while every reader runs
lims <- c(dense = 1e6, sparse = 100)
on <- function(arm, expr) {
  op <- options(frmtmb.joint_root_max = lims[[arm]])
  on.exit(options(op))
  expr
}
fits <- list()
for (arm in names(lims)) {
  fits[[arm]] <- on(arm, {
    f <- quiet(frm(fo, data = d, REML = TRUE))
    invisible(quiet(frmtmb:::sdr_of(f)))
    f
  })
}
fd <- fits$dense
fs <- fits$sparse
cat("joint parameters:", nrow(frmtmb:::sdr_of(fd)$jointPrecision), "\n")
cat("same estimates:", identical(fd$opt$par, fs$opt$par),
    identical(fd$estimates, fs$estimates), "\n")
cat("lost:", paste(frmtmb:::sdr_of(fd)$se_lost), "\n")
cat("root dense:", dim(frmtmb:::sdr_of(fd)$joint_root),
    " sparse:", dim(frmtmb:::sdr_of(fs)$joint_root), "\n")
cat("route dense:", on("dense", class(frmtmb:::get_joint_cov(fd)$V)[1]),
    " sparse:", on("sparse", class(frmtmb:::get_joint_cov(fs)$V)[1]), "\n")
rm("Vjoint", envir = fd$cache)
rm("Vjoint", envir = fs$cache)
rel <- function(x, y) {
  if (!length(x) || !length(y)) return(c(max = NA, scaled = NA))
  x <- as.numeric(as.matrix(x))
  y <- as.numeric(as.matrix(y))
  s <- max(abs(x), abs(y))
  c(max = max(abs(x - y) / pmax(abs(x), abs(y), 1e-300)),
    scaled = max(abs(x - y)) / s)
}
nd <- data.frame(x = seq(1, 2, by = 0.1), g = d$g[1])
for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
nd2 <- nd
nd2$g <- d$g[seq_len(nrow(nd))]
readers <- list(
  fitted_newdata_noRE = function(f) {
    stats::fitted(f, newdata = nd, re_formula = NA)[, "Est.Error"]
  },
  fitted_newdata_RE = function(f) {
    stats::fitted(f, newdata = nd2)[, "Est.Error"]
  },
  fitted_insample = function(f) stats::fitted(f)[, "Est.Error"],
  linpred_se = function(f) {
    frm_linpred(f, newdata = nd2, se.fit = TRUE)$se.fit
  },
  vcov = function(f) stats::vcov(f),
  vcov_full = function(f) stats::vcov(f, full = TRUE),
  draw_space_V = function(f) frmtmb:::fit_draw_space(f)$V,
  draw_space_rootrows = function(f) {
    tcrossprod(frmtmb:::fit_draw_space(f)$root)
  },
  hypothesis = function(f) {
    h <- hypothesis(f, c("x1 + x2 = 0", "x7 = 0", "Intercept = 0"))
    h$hypothesis$Est.Error
  }
)
for (nm in names(readers)) {
  t1 <- system.time(vd <- on("dense", quiet(readers[[nm]](fd))))[["elapsed"]]
  t2 <- system.time(vs <- on("sparse", quiet(readers[[nm]](fs))))[["elapsed"]]
  r <- rel(vd, vs)
  cat(sprintf("%-22s n=%6d  rel max %.3e  scaled %.3e  dense %.2f s  sparse %.2f s\n",
              nm, length(as.numeric(as.matrix(vd))), r[["max"]],
              r[["scaled"]], t1, t2))
}
