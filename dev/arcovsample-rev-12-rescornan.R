# REVIEW script 12, claim 2's downstream consequence: on a Student-t
# set_rescor(TRUE) model with NO autocorrelation (so the REFERENCE build
# accepts it), how many of the draws give a NaN log_lik()?
#
#   Rscript dev/arcovsample-rev-12-rescornan.R <lane|ref>
#
# Own seed (777). The draws are the SAME on both builds because the
# parameter matrix is fixed here, not sampled: nu is walked across the
# range a weakly identified nu reaches.

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]
LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
newx <- "arma_cond_resp" %in% getNamespaceExports("frmtmb")
cat("ARM ", arm, " newexports=", newx, "\n", sep = "")
stopifnot(identical(newx, identical(arm, "lane")))

set.seed(777L)
N <- 50L
x <- rnorm(N)
E <- matrix(rnorm(2 * N), N) %*% chol(matrix(c(1, 0.4, 0.4, 1), 2)) /
  sqrt(rchisq(N, 5) / 5)
dd <- data.frame(x = x, y1 = 0.4 + 0.9 * x + E[, 1],
                 y2 = -0.1 + 0.5 * x + E[, 2])
fit <- frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + student(),
           data = dd, dry_run = "objective")
lab <- c(frmtmb::brms_par_labels(fit), "lp__")
cat("labels: ", paste(lab, collapse = " "), "\n", sep = "")
nus <- c(3, 10, 100, 1e3, 1e5, 1e7, 1e9, 1e11, 1e15, 1e20, 1e50, 1e100,
         1e200, 1e300, 1e306)
nd <- length(nus)
m <- matrix(NA_real_, nd, length(lab), dimnames = list(NULL, lab))
fill <- list(b_y1_Intercept = 0.4, b_y1_x = 0.9, b_y2_Intercept = -0.1,
             b_y2_x = 0.5, sigma_y1 = 1.0, sigma_y2 = 0.9, nu = nus,
             lp__ = 0)
for (nm in names(fill)) m[, nm] <- fill[[nm]]
rc <- grep("^(rescor|thetar)", lab, value = TRUE)
for (nm in rc) m[, nm] <- 0.4
un <- lab[colSums(is.na(m)) > 0L]
stopifnot(length(un) == 0L)
ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                class = "frmtmb_draws")
ll <- log_lik(ds)
rs <- rowSums(ll)
cat("\n", sprintf("%10s %18s %8s", "nu", "rowSum log_lik", "NaN"), "\n",
    sep = "")
for (k in seq_len(nd)) {
  cat(sprintf("%10.1e %18.8f %8d\n", nus[k], rs[k], sum(!is.finite(ll[k, ]))))
}
cat("\ndraws with a non-finite cell: ", sum(!is.finite(rs)), " of ", nd,
    "\n", sep = "")
cat("draws whose rowSum is more than 1 log unit from the nu = 1e300 one: ",
    sum(abs(rs - rs[nd - 1L]) > 1 & nus > 1e7, na.rm = TRUE), " of ",
    sum(nus > 1e7), "\n", sep = "")
lo <- tryCatch(suppressWarnings(loo(ds)), error = identity)
cat("loo(): ", if (inherits(lo, "condition")) conditionMessage(lo) else
  format(lo$estimates["elpd_loo", "Estimate"], digits = 8), "\n",
  sep = "")
cat("DONE\n")
