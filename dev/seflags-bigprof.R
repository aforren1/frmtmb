# Lane seflags: where the SE check spends its time on the 2050-level
# sparse bernoulli fit of dev/seflags-bigp.R (Rprof, by function).
#   Rscript dev/seflags-bigprof.R <lane|base>
a <- commandArgs(trailingOnly = TRUE)
libs <- switch(a[1],
  lane = c("C:/Users/adf44/source/r/wt-seflags-lib",
           "C:/Users/adf44/source/r/rellib-r7"),
  base = "C:/Users/adf44/source/r/rellib-r7")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); ns <- asNamespace("frmtmb")
set.seed(1); L <- 2050
d <- data.frame(f = factor(rep(seq_len(L), each = 4)), x = rnorm(4 * L))
d$y <- rbinom(4 * L, 1, plogis(0.5 + 0.5 * d$x))
f <- suppressWarnings(frm(y ~ f + x, family = bernoulli(), data = d,
  control = frmtmb_control(sparse_x = TRUE, check_se = "ignore")))
pf <- tempfile()
Rprof(pf, interval = 0.05)
invisible(suppressWarnings(ns$se_check(f, frmtmb_control(check_se = "warning"))))
Rprof(NULL)
s <- summaryRprof(pf)
print(head(s$by.total[grepl("se_|cov_|eigen|solve|he|%[*]%|crossprod|colSums|apply", rownames(s$by.total)), ], 30))
