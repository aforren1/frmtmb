# Review of lane gpselimit, check 3: gp(by = f) on sweep rows 3 and 21
# (f alternating a, b), the response moved by y * (1 + 2e-15 * j) for
# j = 0..20: how many white sub-GP hyperparameters keep a finite SE, the
# lost reasons, the conditions, and the SEs of the intercept and sigma.
# Usage: Rscript gpselimit-review-byjit.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-review-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")
for (row in c(3, 21)) for (j in 0:20) {
  s <- rv_sweep_row(row, j)
  s$d$f <- factor(rep(c("a", "b"), length.out = nrow(s$d)))
  q <- rv_fit(frm(y ~ gp(x, by = f, k = s$k), data = s$d))
  fit <- q$fit
  ss <- rv_ses(fit)
  thp <- which(names(ss$se) == "theta")
  wh <- integer(0)
  for (bk in rv_hsgp(fit)) {
    if (ns$gp_basis_white(bk, fit$estimates$theta[bk$theta_idx])) {
      wh <- c(wh, thp[bk$theta_idx])
    }
  }
  cat(sprintf(paste0("BYJIT %s row %d j %2d | white theta %d finite SE %d | ",
                     "beta %.10g sigma %.10g | lost [%s] | conds [%s]\n"),
              arm, row, j, length(wh), sum(is.finite(ss$se[wh])),
              ss$se[["beta"]], ss$se[["betad"]],
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              rv_cstr(q$conds)))
}
