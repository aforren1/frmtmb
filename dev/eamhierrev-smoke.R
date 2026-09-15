# REVIEW of lane eamhier: verify the restored reference build by FITTING
# rather than by reading a version string, per the round's instruction
# that a hollow directory still answers packageVersion().
#
# IT DOES NOT REPRODUCE dev/machine-library.md's two restore numbers,
# and that is a fact about that file rather than about the library:
# it records -268.5330674 and -128.3133341 without the script that
# made them, so the data-generating code below is a GUESS at the
# construction and produces -302.678 and -215.669 instead. Recording
# the construction beside those two numbers would close this.
#
# The instrument check that DID settle the library is in
# dev/eamhierrev-sv.R: single-level replicate seed 20260910 refits to
# this lane's own recorded logLik of -6657.979468 in every digit.
#
# Run: Rscript --vanilla dev/eamhierrev-smoke.R

.libPaths(c("C:/Users/adf44/source/r/eamhierrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("frmtmb ", format(packageVersion("frmtmb")),
    " from ", dirname(find.package("frmtmb")), "\n", sep = "")
cat("frmtmb.eam ", format(packageVersion("frmtmb.eam")),
    " from ", dirname(find.package("frmtmb.eam")), "\n", sep = "")
cat("StanHeaders ", format(packageVersion("StanHeaders")), "\n",
    sep = "")

set.seed(1)
n <- 200L
g <- factor(rep(seq_len(20L), each = 10L))
u <- stats::rnorm(20L, 0, 0.7)
x <- stats::rnorm(n)
y <- 1 + 0.5 * x + u[g] + stats::rnorm(n, 0, 1)
d <- data.frame(y = y, x = x, g = g)
f1 <- frm(y ~ x + (1 | g), family = gaussian(), data = d, se = FALSE)
cat(sprintf("gaussian GLMM logLik %.7f  (reference -268.5330674)\n",
            as.numeric(stats::logLik(f1))))

set.seed(2)
dd <- ddm_simulate(400L, mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.5)
f2 <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = wiener(), data = dd, se = FALSE)
cat(sprintf("wiener DDM logLik  %.7f  (reference -128.3133341)\n",
            as.numeric(stats::logLik(f2))))
