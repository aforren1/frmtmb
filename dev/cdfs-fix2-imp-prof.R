# Fix round 2, N4: where the importance-sampled zero-truncated poisson
# spends its time, by Rprof, on the given library.
#   Rscript dev/cdfs-fix2-imp-prof.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", find.package("frmtmb"), "\n")
set.seed(21)
ng <- 40L
per <- 6L
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(ng * per)
u <- rnorm(ng, 0, 0.8)
y <- rpois(ng * per, exp(1.1 + 0.4 * x + u[g]))
dd <- data.frame(y = y, x = x, g = g)[y >= 1L, ]
bfm <- bf(y | trunc(lb = 1) ~ x + (1 | g)) + poisson()
pf <- tempfile()
Rprof(pf, interval = 0.01)
fi <- frm(bfm, dd, importance = 500L)
Rprof(NULL)
s <- summaryRprof(pf)
cat("== by self time\n")
print(utils::head(s$by.self, 15))
cat("== by total time (frmtmb functions)\n")
bt <- s$by.total
print(utils::head(bt[grepl("row_lpdf|log_window|fam_log_pair|poisson_cdf|log1m|trunc_edges|cdf_edges|MakeADFun|imp_|importance|optim|nlminb", rownames(bt)), ], 25))
cat("END\n")
