# Short check of lane sparsese: on a lost-sd fit (100 groups, group sd 0,
# seed 1, the indep design), frm_lp_basis()'s V_root V_root' + V_rest on
# the dense route (limit 1e6) against the sparse route (limit 20), and the
# size of the b-row variances against the rows' variances.
#   Rscript sparsese-short-lostb.R <seed>
seed <- as.integer(commandArgs(trailingOnly = TRUE)[1])
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "seed", seed, "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
ng <- 100L
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, 0)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
nd2 <- data.frame(x = seq(1, 2, by = 0.1), g = d$g[1:11])
for (k in 1:7) nd2[[paste0("x", k)]] <- nd2$x^k
R <- list()
for (lim in c(1e6, 20)) {
  options(frmtmb.joint_root_max = lim)
  fit <- q(frm(fo, data = d, REML = TRUE))
  s <- q(frmtmb:::sdr_of(fit))
  lb <- q(frm_lp_basis(fit, newdata = nd2))
  cat("limit", lim, "lost:", paste(names(s$se_lost), s$se_lost),
      " theta", format(fit$opt$par[["theta"]], digits = 6),
      " route", class(q(frmtmb:::get_joint_cov(fit))$V)[1], "\n")
  rn <- rownames(s$jointPrecision)[lb$coef_pos]
  Rm <- tcrossprod(lb$V_root) + lb$V_rest
  dimnames(Rm) <- list(rn, rn)
  A <- as.matrix(lb$A)
  R[[as.character(lim)]] <- list(R = Rm, A = A, cp = lb$coef_pos,
                                 v = rowSums((A %*% lb$V_root)^2))
}
a <- R[[1]]; b <- R[[2]]
# the dense route keeps every column; the sparse one only those loaded
keep <- match(b$cp, a$cp)
a$R <- a$R[keep, keep]
cat("positions: dense", length(a$cp), " sparse", length(b$cp), "
")
ib <- which(rownames(a$R) == "b")
s <- sqrt(diag(a$R))
cat(sprintf("dense vs sparse V_root V_root' + V_rest: scaled max %.3e;",
            max(abs(a$R - b$R) / outer(s, s))),
    sprintf("b diag rel max %.3e\n",
            max(abs(diag(a$R)[ib] - diag(b$R)[ib]) / diag(a$R)[ib])))
cat("b-row variances (dense):", format(range(diag(a$R)[ib]), digits = 3),
    "; rows' variances:", format(range(a$v), digits = 3), "\n")
cat(sprintf("row variances dense vs sparse: max rel %.3e\n",
            max(abs(a$v - b$v) / a$v)))
