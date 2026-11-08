# Release review 0.80.0: which start reaches the maximum (-638.885) of the
# binomial mm 22 row? A grid of starts, se = FALSE.
#   Rscript dev/review-rel-mm22-starts.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
set.seed(22001)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                xq = round(runif(n, 0, 5), 3),
                xb = round(runif(n, 0, 3000), 1),
                g = factor(sample(1:15, n, TRUE)))
d$ug <- rnorm(15, 0, 0.4)[d$g]
d$nt <- 10L
d$y <- rbinom(n, 10, plogis(3 * d$xp / (1 + d$xp) - 2))
fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1, K ~ 1, nl = TRUE)
for (vm in as.numeric(strsplit(Sys.getenv("VMS", "-3,-1,-0.3,0.3,1,3"), ",")[[1]])) for (k in as.numeric(strsplit(Sys.getenv("KS", "-3,-0.5,-0.1,0.5,3"), ",")[[1]])) {
  f <- tryCatch(suppressWarnings(suppressMessages(frm(fm, family = binomial(),
    data = d, se = FALSE, start = list(beta = c("Vm_(Intercept)" = vm,
                                                "K_(Intercept)" = k))))),
    error = function(e) NULL)
  if (is.null(f)) next
  cat(sprintf("start Vm %5.1f K %5.1f -> ll %.4f Vm %.5g K %.5g\n", vm, k,
              as.numeric(logLik(f)), f$estimates$beta[[1]],
              f$estimates$beta[[2]]))
}
cat("STARTS DONE\n")
