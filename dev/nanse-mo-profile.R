# Defect 8: is frmtmb's optimum the maximum likelihood? Given the two
# simplexes, `ls ~ mo(income) * age` is a linear model, so the profile
# log-likelihood over the simplexes is exact (OLS, sigma^2 = RSS / n).
# Maximized over stick-breaking coordinates in [0, 1]^4 with box bounds,
# which put the simplex boundary at a finite point, from 9 starts.
# Compared with frmtmb's logLik per seed.
#   Rscript dev/nanse-mo-profile.R [lib] [seeds] [out]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:200
out <- if (length(args) > 2) args[3] else "dev/nanse-log/mo-profile.tsv"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
stick <- function(s) c(1 - s[1], s[1] * (1 - s[2]), s[1] * s[2])
prof <- function(s, d) {
  w1 <- stick(s[1:2])
  w2 <- stick(s[3:4])
  code <- as.integer(d$income) - 1L
  # brms's mo(): D * sum of the first `code` simplex weights, D = 3
  m1 <- 3 * c(0, cumsum(w1))[code + 1L]
  m2 <- 3 * c(0, cumsum(w2))[code + 1L]
  X <- cbind(1, d$age, m1, m2 * d$age)
  r <- stats::lm.fit(X, d$ls)$residuals
  n <- nrow(d)
  -(n / 2) * (log(2 * pi * sum(r^2) / n) + 1)
}
rows <- list()
st <- as.matrix(expand.grid(a = c(0.2, 0.5, 0.8), b = c(0.2, 0.8)))
for (s in seeds) {
  d <- mk(s)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  best <- -Inf
  bpar <- NULL
  for (i in seq_len(nrow(st))) {
    for (j in c(0.3, 0.7)) {
      p0 <- c(st[i, 1], st[i, 2], j, 1 - j)
      o <- stats::nlminb(p0, function(p) -prof(p, d), lower = 0, upper = 1)
      if (-o$objective > best) {
        best <- -o$objective
        bpar <- o$par
      }
    }
  }
  w2 <- stick(bpar[3:4])
  w1 <- stick(bpar[1:2])
  zt <- f$estimates[grepl("^zeta", names(f$estimates))]
  fw <- lapply(zt, function(z) {
    x <- exp(c(0, z))
    x / sum(x)
  })
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, ll_frmtmb = as.numeric(logLik(f)), ll_profile = best,
    gap = best - as.numeric(logLik(f)),
    w1_prof = paste(signif(w1, 3), collapse = ","),
    w2_prof = paste(signif(w2, 3), collapse = ","),
    w1_frm = paste(signif(fw[[1]], 3), collapse = ","),
    w2_frm = paste(signif(fw[[2]], 3), collapse = ","),
    prof_boundary = sum(c(w1, w2) < 1e-6))
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("seeds", nrow(X), "\n")
cat("gap = profile max - frmtmb logLik: quantiles\n")
print(stats::quantile(X$gap, c(0, 0.5, 0.9, 0.99, 1)))
cat("seeds with gap > 1e-3:", sum(X$gap > 1e-3), "; > 1e-2:",
    sum(X$gap > 1e-2), "; > 0.1:", sum(X$gap > 0.1), "\n")
cat("profile optimum on the simplex boundary (a weight < 1e-6):",
    sum(X$prof_boundary > 0), "\n")
