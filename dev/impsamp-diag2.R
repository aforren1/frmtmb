# Diagnosis 2: why does the re-anchored iteration diverge on p124 with
# ESS 0.998? Candidate: the frozen tape's CURVATURE along directions the
# conditional mode must follow is a noisy sample statistic. Measure, at
# the Laplace optimum on the base build:
#  (a) weight variance, ESS and MCSE per group at the anchor;
#  (b) the frozen tape's gradient against a central difference of the
#      honest value (the re-anchored objective, quadrature where exact);
#  (c) the frozen tape's Hessian against the exact (quadrature) Hessian,
#      over 20 draw seeds, and the contraction factor of the fixed-point
#      map, the spectral radius of I - Hn^-1 H_exact;
#  (d) the curvature's prediction from theory: for intercept beta,
#      sum_g I_g (1 - mean z_g^2) + marginal information.
# Seeds: data 124; draw seeds 1:20. Log dev/impsamp-out/diag2.txt
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")
E <- new.env(parent = ns)
sys.source("dev/equiv-specs.R", envir = E)
db <- E$eq_dat("binomial", seed = 124)
form <- bf(y | trials(ntr) ~ x + (1 | g))
lap <- frm(form, family = binomial(), data = db)
quad <- frm(form, family = binomial(), data = db, quadrature = TRUE)
fr <- lap$frame
tpl <- lap$obj$env$parList(lap$opt$par)
nll <- ns$build_objective(fr)
lobj <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
lay <- ns$imp_layout(fr[["re_blocks"]])
gm <- ns$imp_group_map(fr, lay)
otpl <- ns$imp_template(tpl, "b")
freeze <- function(at, N, seed) {
  plan <- ns$imp_plan(lobj, fr, lay, at, N, seed)
  io <- ns$build_importance_objective(fr, lay, gm, plan)
  pl <- if (length(formals(ns$imp_par_list)) >= 3L) ns$imp_par_list(tpl, at, fr[["map"]]) else ns$imp_par_list(tpl, at)
  am <- as.matrix(io$amat(pl))
  io$set_offset(apply(am, 1L, max))
  list(plan = plan, io = io, am = am,
       obj = RTMB::MakeADFun(io$fn, otpl, silent = TRUE))
}
at <- lap$opt$par
cat("Laplace par:", format(at, digits = 8), "\n")
cat("Quadrature par:", format(quad$opt$par, digits = 8), "\n")
numh <- function(gr, p, h = 1e-4) {
  m <- sapply(seq_along(p), function(j) {
    e <- replace(numeric(length(p)), j, h)
    (gr(p + e) - gr(p - e)) / (2 * h)
  })
  0.5 * (m + t(m))
}
Hq <- numh(quad$obj$gr, at)
Hl <- numh(lobj$gr, at)
cat("H exact (quadrature tape) at the Laplace optimum:\n")
print(signif(Hq, 5))
cat("H Laplace at the Laplace optimum:\n")
print(signif(Hl, 5))
# conditional information of each group, from the Hessian's diagonal
full <- lobj$env$last.par
lobj$fn(at)
full <- lobj$env$last.par
hd <- Matrix::diag(lobj$env$spHess(full, random = TRUE))
sig2 <- exp(2 * at[["theta"]])
cat("per-group conditional precision H_g:", format(hd, digits = 4), "\n")
cat("prior precision 1 / sigma^2:", format(1 / sig2, digits = 4),
    " H_g sigma^2 mean:", format(mean(hd) * sig2, digits = 4), "\n")
for (N in c(200L, 512L, 1000L, 4000L)) {
  cat("\n== N", N, "==\n")
  rows <- list()
  for (s in 1:20) {
    fz <- freeze(at, N, s)
    es <- ns$imp_ess(fz$am, N)
    w <- exp(fz$am - apply(fz$am, 1L, max))
    Hn <- fz$obj$he(at)
    g <- fz$obj$gr(at)
    gq <- quad$obj$gr(at)
    ev <- eigen(solve(Hn, Hq), only.values = TRUE)$values
    rho <- max(Mod(1 - ev))
    # theory: the (beta_0, beta_0) entry
    z <- ns$imp_draws(lay[["nb"]], N, s)
    mz2 <- rowMeans(z^2)
    th <- sum((hd - 1 / sig2) * (1 - mz2)) + Hq[1, 1]
    rows[[s]] <- data.frame(seed = s, ess_min = min(es$ess), mcse = es$mcse,
                            sdlogw = max(apply(fz$am, 1L, stats::sd)),
                            gmax_err = max(abs(g - gq)),
                            H11 = Hn[1, 1], H11_theory = th,
                            mineig = min(eigen(Hn, TRUE, TRUE)$values),
                            rho = rho)
  }
  r <- do.call(rbind, rows)
  print(format(r, digits = 4), row.names = FALSE)
  cat(sprintf("rho > 1 in %d of 20 seeds; median rho %.3f; H11 exact %.3f;",
              sum(r$rho > 1), median(r$rho), Hq[1, 1]),
      sprintf("sd(H11) %.3f; cor(H11, theory) %.4f\n", sd(r$H11),
              cor(r$H11, r$H11_theory)))
}
