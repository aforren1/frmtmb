# Memory of gp_krig_cov() under gp(x, by = <numeric>), on a 2000-row
# grid (the review's m9, dev/cifixrev-krigmem.R measured peak RSS and
# could not tell it from load). Two load-independent instruments:
# Rprofmem's large allocations inside the call (each n x n is 32 MB),
# and gc()'s "max used" Vcells after a reset, in units of n^2 doubles.
# Usage: Rscript dev/ciharden-krigmem.R <base|lane>
arm <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/wt-ciharden-lib"
.libPaths(c(if (arm == "lane") lib, "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(5)
n0 <- 60
d <- data.frame(x = round(stats::runif(n0, 0, 6), 1),
                w = stats::runif(n0, 0.5, 2))
d$y <- 0.5 + sin(d$x) * d$w + stats::rnorm(n0, 0, 0.3)
fw <- frm(bf(y ~ gp(x, by = w)), data = d)
m <- 2000
nd <- data.frame(x = seq(-0.5, 6.5, length.out = m) + 1e-3,
                 w = stats::runif(m, 0.5, 2))
ed <- frmtmb:::lp_eta_design(fw, fw$frame$linpreds[["y.mu"]], nd, FALSE,
                             FALSE)
sp <- Filter(function(s) !is.null(s$krig), ed$sm_parts)[[1L]]
kg <- sp$krig
f <- frmtmb:::gp_krig_cov
cat("profmem capability:", capabilities("profmem"), "\n")
tf <- tempfile()
invisible(gc())
Rprofmem(tf, threshold = m^2 * 8 / 2)
S <- f(kg)
Rprofmem(NULL)
lns <- grep("^[0-9]", readLines(tf), value = TRUE)
print(substr(lns, 1, 200))
big <- as.numeric(sub(" *:.*", "", grep("^[0-9]", readLines(tf),
                                          value = TRUE)))
cat(sprintf("allocations of at least n^2/2 doubles: %d, total %.2f n^2\n",
            length(big), sum(big) / (8 * m^2)))
rm(S)
invisible(gc())
base_vc <- gc(reset = TRUE)["Vcells", "used"]
S <- f(kg)
mx <- gc()["Vcells", "max used"]
cat(sprintf("peak Vcells above the start: %.2f n^2\n",
            (mx - base_vc) / m^2))
