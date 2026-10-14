# Which build a re-measure ran against, verified by a fit and not by a
# version string (a hollow package directory still answers
# packageVersion()). Run by dev/remeasure-run.sh with PORT_LIB set.
#
#   PORT_LIB=<libs> Rscript dev/remeasure-provenance.R
lib <- Sys.getenv("PORT_LIB", unset = "")
libs <- c(strsplit(lib, ";", fixed = TRUE)[[1]],
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
cat("libPaths:", paste(.libPaths(), collapse = " ; "), "\n")
for (l in strsplit(lib, ";", fixed = TRUE)[[1]]) {
  d <- list.dirs(l, recursive = FALSE, full.names = FALSE)
  d <- setdiff(d, grep("^00LOCK", d, value = TRUE))
  hollow <- d[!file.exists(file.path(l, d, "DESCRIPTION"))]
  cat(sprintf("%s: %d package dirs, %d hollow %s\n", l, length(d),
              length(hollow), paste(hollow, collapse = " ")))
}
for (p in c("frmtmb", "frmtmb.sample", "frmtmb.eam", "frmtmb.latent",
            "frmtmb.spline", "frmtmb.coupling", "frmtmb.learn",
            "frmtmb.ode", "brms", "rstan", "tmbstan", "RTMB")) {
  where <- tryCatch(dirname(find.package(p)), error = function(e) NA)
  cat(sprintf("%-16s %-8s %s\n", p,
              if (is.na(where)) "-" else format(packageVersion(p, where)),
              if (is.na(where)) "not installed" else where))
}
suppressMessages(library(frmtmb))
# The fit dev/correct-libcheck.R records (seed 1, n = 200, 10 groups):
# gaussian -295.602189818 and poisson -332.137876329 on the restored
# library of 2026-09-28.
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 20)))
d$y <- 1 + 0.5 * d$x + rnorm(10, 0, 0.7)[d$g] + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(10, 0, 0.5)[d$g]))
f1 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
f2 <- frm(bf(cnt ~ x + (1 | g)), family = poisson(), data = d)
cat(sprintf("fit check: gaussian (1 | g) logLik %.9f, poisson %.9f\n",
            as.numeric(logLik(f1)), as.numeric(logLik(f2))))
