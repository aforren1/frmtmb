# Lane squant, fix round 1: the frm_compat() rows of quantile_elf against
# asym_laplace's, with every extension loaded, and the fits behind the
# rows: trunc() with fitted(), cens(), mi() on the response and on a
# predictor, a mixture of two ELFs.
# Run: Rscript dev/squant-compat.R > dev/squant-out/compat.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  for (p in c("frmtmb.sample", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
              "frmtmb.ode", "frmtmb.survival", "frmtmb.spatial",
              "frmtmb.coupling")) library(p, character.only = TRUE)
})
cat("lib:", find.package("frmtmb"), "\n")
a <- frm_compat("asym_laplace")
q <- frm_compat("quantile_elf")
key <- function(x) paste(x$feature_a, x$feature_b)
names(x <- q)
m <- merge(data.frame(k = key(a), al = a$status),
           data.frame(k = key(q), elf = q$status), by = "k", all = TRUE)
m$k <- sub("asym_laplace|quantile_elf", "FAM", m$k)
mm <- stats::aggregate(cbind(al, elf) ~ k, m, function(v) v[1],
                       na.action = na.pass)
cat("pairs: asym_laplace", nrow(a), "quantile_elf", nrow(q), "\n")
d <- mm[is.na(mm$al) | is.na(mm$elf) | mm$al != mm$elf, ]
print(d, row.names = FALSE)
for (f in c("trunc()", "cens()", "mi()", "mi_pred()", "mixture",
            "frm_sample", "quadrature", "REML")) {
  r <- frm_compat("quantile_elf", f)
  cat(f, ":", r$status, "|", substr(r$note, 1, 150), "\n")
}
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 1 + d$x + rnorm(n)
d$yp <- abs(d$y) + 0.1
tf <- frm(bf(yp | trunc(lb = 0) ~ x, quantile = 0.7),
          family = quantile_elf(), data = d)
cat("trunc fit conv", tf$opt$convergence, "\n")
print(tryCatch(head(fitted(tf), 2), error = function(e) conditionMessage(e)))
d$ym <- d$y
d$ym[1:10] <- NA
mf <- frm(bf(ym | mi() ~ x, quantile = 0.5), family = quantile_elf(),
          data = d)
cat("mi response conv", mf$opt$convergence, "\n")
d$xm <- d$x
d$xm[1:10] <- NA
pf <- frm(bf(y ~ mi(xm), quantile = 0.5) +
            bf(xm | mi() ~ 1, quantile = 0.5) + set_rescor(FALSE),
          family = quantile_elf(), data = d)
cat("mi predictor conv", pf$opt$convergence, "\n")
