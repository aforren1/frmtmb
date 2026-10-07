# Two unseen gp() positions one bit apart: one kriged position or two?
# The fit of test-gp-multidim.R (seed 17, 120 rows on a 0.1 grid), new
# rows u, u * (1 + 2^-52), u with u = 10/3. Prints the kriging
# covariance of the three rows and the gap (S11 - S12) in units of the
# nugget's share, gp_nugget * sd^2: 0 when the two are merged, 1 when
# they are two positions.
# Usage: Rscript dev/ciharden-poskey.R <lib or base>
a <- commandArgs(TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(17)
xg <- round(runif(120, 0, 10), 1)
dg <- data.frame(y = sin(xg) + rnorm(120, 0, 0.3), x = xg)
fg <- frm(bf(y ~ gp(x)) + gaussian(), data = dg)
u <- 10 / 3
nd <- data.frame(x = c(u, u * (1 + 2^-52), u))
ed <- frmtmb:::lp_eta_design(fg, fg$frame$linpreds[["y.mu"]], nd, FALSE,
                             FALSE)
kg <- Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1L]]
S <- frmtmb:::gp_krig_cov(kg)
print(S, digits = 17)
cat(sprintf("gap (S11 - S12) / (gp_nugget sd^2) = %.12g\n",
            (S[1, 1] - S[1, 2]) / (frmtmb:::gp_nugget * kg$sd2)))
lb <- frm_lp_basis(fg, newdata = nd)
A <- unname(as.matrix(lb$A))
cat("design rows 1 and 2 identical:", identical(A[1, ], A[2, ]),
    "; max |difference|:", max(abs(A[1, ] - A[2, ])), "\n")
