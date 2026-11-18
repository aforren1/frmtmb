# rel084-kt copy of dev/nbrobust-jac.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Lane nbrobust: value and derivative in the log shape of the negbinomial
# log mass, on the tape, against a 256-bit reference.
#   Rscript dev/nbrobust-jac.R [out.tsv]
# Arms: "robust" is RTMB::dnbinom_robust(y, lmu, 2 lmu - lsh); "scalar"
# and "vector" are frmtmb:::nb_lpmf(y, lmu, lsh) with one shape for all
# rows and with a shape per row (the two branches of nb_rising()).
# Grid: log shape 0 to 40 by 0.5, counts 0 to 1e4, means 0.1 to 1e4.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "dev/nbrobust-out/jac.tsv"
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB)
  library(Rmpfr)
})
cat("lib:", dirname(system.file(package = "frmtmb")), " version",
    as.character(packageVersion("frmtmb")), " RTMB",
    as.character(packageVersion("RTMB")), " TMB",
    as.character(packageVersion("TMB")), "\n")
nb_lpmf <- frmtmb:::nb_lpmf
B <- 256
ys <- c(0, 1, 2, 3, 5, 10, 24, 25, 30, 100, 1000, 1e4)
mus <- c(0.1, 3, 100, 1e4)
lshs <- seq(0, 40, by = 0.5)
ref <- function(y, lmu, lsh) {
  n <- exp(mpfr(lsh, B))
  mu <- exp(mpfr(lmu, B))
  v <- lgamma(y + n) - lgamma(n) - lgamma(mpfr(y + 1, B)) +
    n * log(n / (n + mu)) + y * log(mu / (n + mu))
  dg <- if (y == 0) mpfr(0, B) else digamma(y + n) - digamma(n)
  dsh <- n * (dg - log1p(mu / n) + (mu - y) / (n + mu))
  dmu <- y - (n + y) * mu / (n + mu)
  c(v = asNumeric(v), dsh = asNumeric(dsh), dmu = asNumeric(dmu))
}
rows <- list()
for (mu in mus) {
  lmu <- log(mu)
  N <- length(ys)
  tp <- list(
    robust = MakeTape(function(p) {
      RTMB::dnbinom_robust(ys, p[1] + 0 * ys, 2 * p[1] - p[2], log = TRUE)
    }, c(lmu, 1)),
    scalar = MakeTape(function(p) {
      nb_lpmf(ys, p[1] + 0 * ys, p[2])
    }, c(lmu, 1)),
    vector = MakeTape(function(p) {
      nb_lpmf(ys, p[1] + 0 * ys, p[2] + 0 * ys)
    }, c(lmu, 1))
  )
  jf <- lapply(tp, function(t) t$jacfun())
  for (lsh in lshs) {
    rf <- vapply(ys, function(y) ref(y, lmu, lsh), numeric(3))
    for (arm in names(tp)) {
      v <- tp[[arm]](c(lmu, lsh))
      J <- jf[[arm]](c(lmu, lsh))
      rows[[length(rows) + 1L]] <- data.frame(
        arm = arm, mu = mu, lsh = lsh, y = ys, v = v, v_ref = rf["v", ],
        dsh = J[, 2], dsh_ref = rf["dsh", ], dmu = J[, 1],
        dmu_ref = rf["dmu", ])
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, out, sep = "\t", row.names = FALSE, quote = FALSE)
eps <- .Machine$double.eps
tab$sc <- eps * (tab$y + tab$mu + 1)
tab$e_dsh <- abs(tab$dsh - tab$dsh_ref)
tab$e_dmu <- abs(tab$dmu - tab$dmu_ref)
tab$e_v <- abs(tab$v - tab$v_ref) / pmax(abs(tab$v_ref), 1e-300)
band <- cut(tab$lsh, c(-1, 10, 18, 24, 30, 35, 40),
            labels = c("0-10", "10.5-18", "18.5-24", "24.5-30", "30.5-35",
                       "35.5-40"))
cat("\nmax |error| of the derivative in the log shape, by log-shape band\n")
cat("and in units of eps * (y + mu + 1); the value's relative error\n")
for (arm in c("robust", "scalar", "vector")) {
  for (b in levels(band)) {
    s <- tab[tab$arm == arm & band == b, ]
    cat(sprintf(paste("%-6s lsh %-8s n %4d  max|d dsh| %.3e  max units",
                      "%9.3g  max|d dmu| %.3e  max rel value %.3e\n"),
                arm, b, nrow(s), max(s$e_dsh), max(s$e_dsh / s$sc),
                max(s$e_dmu), max(s$e_v)))
  }
}
# the value against dnbinom_robust where that one is right (relative
# error to the reference below 1e-13)
r <- tab[tab$arm == "robust", ]
s <- tab[tab$arm == "scalar", ]
ok <- r$e_v < 1e-13
cat(sprintf(paste("\nnb_lpmf (scalar) against dnbinom_robust where it is",
                  "right: %d of %d points, max rel %.3e\n"),
            sum(ok), nrow(r), max(abs(s$v[ok] - r$v[ok]) /
                                     abs(r$v[ok]))))
cat(sprintf("scalar against vector: max |value| diff %.3e, max |dsh| diff %.3e\n",
            max(abs(s$v - tab$v[tab$arm == "vector"])),
            max(abs(s$dsh - tab$dsh[tab$arm == "vector"]))))
# the full table of the derivative at the item's mean, 3, count 3
cat("\nlog shape, y = 3, mu = 3: robust, nb_lpmf, reference (d / d log shape)\n")
for (l in c(0, 5, 10, 14, 18, 20, 22, 24, 26, 28, 30, 32, 35, 38, 40)) {
  i <- function(a) tab$arm == a & tab$mu == 3 & tab$y == 3 & tab$lsh == l
  cat(sprintf("%5.1f % .5e % .5e % .5e\n", l, tab$dsh[i("robust")],
              tab$dsh[i("scalar")], tab$dsh_ref[i("scalar")]))
}
cat("END\n")
