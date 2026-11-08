# Lane latperf: where fitted() spends its time, and what the fit keeps.
# Usage: Rscript dev/latperf-prof-two.R <arm> <n> <design> [dense]
#   arm and design as in dev/latperf-prof.R. Seed 1, the same generator.
args <- commandArgs(TRUE)
arm <- args[1]
n <- as.integer(args[2])
design <- args[3]
base <- c("C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-latperf-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
if ("dense" %in% args) {
  options(frmtmb.latres_se_rows = Inf, frmtmb.latres_dense_rows = Inf)
}
cat(sprintf("lib: %s frmtmb %s arm %s n %d design %s\n",
            dirname(dirname(system.file(package = "frmtmb"))),
            format(packageVersion("frmtmb")), arm, n, design))
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n),
                id = factor(sample.int(50, n, TRUE)))
e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
u <- rnorm(50, 0, 0.3)
d$y <- rpois(n, exp(.5 + .3 * d$x + e +
                     if (design == "two") u[d$id] else 0))
f <- if (design == "two") {
  y ~ x + (1 | id) + ar(t, g, cov = TRUE)
} else {
  y ~ x + ar(t, g, cov = TRUE)
}
fit <- suppressWarnings(frm(f, family = poisson(), data = d))
tmp <- tempfile()
Rprof(tmp, interval = 0.01)
v <- suppressWarnings(fitted(fit))
Rprof(NULL)
s <- summaryRprof(tmp)
cat(sprintf("sampled %.2f s\n", s$sampling.time))
print(head(s$by.self, 12))
bt <- s$by.total
print(head(bt[order(-bt$total.time), c("total.time", "self.time")], 45))
# what the cache keeps, by entry
for (nm in ls(fit$cache)) {
  cat(sprintf("CACHE %-14s %8.2f MB\n", nm,
              length(serialize(get(nm, envir = fit$cache), NULL)) / 2^20))
}
jc <- fit$cache$Vjoint
if (inherits(jc$V, "frm_joint_sparse")) {
  for (nm in c("L", "J", "Vf", "memo")) {
    cat(sprintf("VJOINT %-6s %8.2f MB\n", nm,
                length(serialize(jc$V[[nm]], NULL)) / 2^20))
  }
}
