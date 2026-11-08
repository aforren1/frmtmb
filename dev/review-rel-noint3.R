# Short check of the 0.80.0 release fixes, item 2: the designs that do not
# span the constant, against lm() (gaussian ML: lm's logLik with REML =
# FALSE), and design_center() returning them unchanged.
#   Rscript dev/review-rel-noint3.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
dc <- get("design_center", asNamespace("frmtmb"))
for (seed in 1:3) for (m in c(1e3, 1e5, 1e6)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(z1 = rnorm(n), z2 = rnorm(n),
                  f = factor(sample(letters[1:3], n, TRUE)))
  d$y <- rnorm(n)
  d$w1 <- d$z1 + m
  d$w2 <- d$z2 + m
  d$y0 <- 0.5 * d$w1 - 0.5 * d$w2 + 0.6 * d$z1 * d$z2 + rnorm(n, 0, 0.5)
  for (fo in list(y0 ~ 0 + w1 + w2, y0 ~ 0 + w1 * w2)) {
    X <- model.matrix(fo, d)
    same <- identical(dc(X), X)
    f <- suppressWarnings(frm(fo, data = d))
    l <- lm(fo, d)
    cat(sprintf("seed %d m %.0e %-18s design_center unchanged %s | frm %.5f lm %.5f diff %.2e lm rank %d\n",
                seed, m, deparse(fo[[3]]), same, as.numeric(logLik(f)),
                as.numeric(logLik(l)),
                as.numeric(logLik(f)) - as.numeric(logLik(l)), l$rank))
  }
}
cat("NOINT3 DONE\n")
