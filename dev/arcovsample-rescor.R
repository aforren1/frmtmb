# Lane wt-arcovsample: log_lik() on a set_rescor(TRUE) model that also
# carries a cov = FALSE ARMA term on each response. The objective shifts
# mu and THEN takes the joint density, so draws_row_loglik() has to do
# the same; rescor_row_loglik() reads the shifted mu because
# arma_cond_dpars() runs before it.
#
#   Rscript dev/arcovsample-rescor.R > dev/arcovsample-log/rescor.txt

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

SEED <- 4023L
set.seed(SEED)
ng <- 8L
nt <- 9L
dd <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
dd <- dd[-c(4L, 20L, 33L), ]
dd$x <- stats::rnorm(nrow(dd))
dd$y <- 1 + 0.5 * dd$x + stats::rnorm(nrow(dd))
dd$y2 <- 0.4 + 0.8 * dd$y + stats::rnorm(nrow(dd), 0, 0.7)
dd <- dd[sample(nrow(dd)), ]
rownames(dd) <- NULL
cat("N =", nrow(dd), " seed =", SEED, "\n")

for (fam in c("gaussian", "student")) {
  f <- if (identical(fam, "gaussian")) gaussian() else student()
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + ar(t, g)) + bf(y2 ~ x + ar(t, g)) +
                 set_rescor(TRUE), family = f, data = dd,
               chains = 1, iter = 400, refresh = 0, seed = 9)))
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  nd <- nrow(ds$draws)
  res <- numeric(nd)
  for (d in seq_len(nd)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, d, idx)
    nll <- as.numeric(frmtmb::build_objective(sh$frame)(sh$estimates))
    res[d] <- sum(ll[d, ]) + nll
  }
  cat("\n---- rescor + ar(1) on both, ", fam, " ----\n", sep = "")
  cat("  draws x cols : ", nd, " x ", ncol(ll), "\n", sep = "")
  cat("  max |row sum + nll| : ", format(max(abs(res)), digits = 12),
      "\n", sep = "")
  cat("  max relative        : ",
      format(max(abs(res) / abs(rowSums(ll))), digits = 12), "\n",
      sep = "")
  # the shift must be in there: the unshifted joint density is a
  # different number by far more than the residual
  sh <- frmtmb.sample:::draws_fit_at(ds, 1L, idx)
  dp0 <- frmtmb::with_cs_offsets(sh, NULL, frmtmb::eval_dpars(sh))
  plain <- sum(frmtmb::rescor_row_loglik(sh, dp0))
  cat("  shifted vs unshifted joint at draw 1 : ",
      format(sum(ll[1, ]), digits = 10), " vs ",
      format(plain, digits = 10), "\n", sep = "")
  lo <- suppressWarnings(loo(ds))
  cat("  elpd_loo : ", format(lo$estimates["elpd_loo", "Estimate"],
                              digits = 10), " (finite ",
      is.finite(lo$estimates["elpd_loo", "Estimate"]), ")\n", sep = "")
}
cat("\nDONE\n")
