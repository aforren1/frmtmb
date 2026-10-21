.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
E <- new.env(parent = asNamespace("frmtmb")); sys.source("R/families-cdf.R", envir = E)
a <- -pi; b <- 1; k <- 1
ks <- get("ad_min", E)(k, 10.5); print(c(E$vm_series_H(b, ks), E$vm_series_H(a, ks)))
kh <- get("ad_max", E)(k, 9); wa <- E$vm_wrap(a); wb <- E$vm_wrap(b); print(unlist(wa)); print(unlist(wb))
za <- E$vm_hill_z(wa$t, kh); zb <- E$vm_hill_z(wb$t, kh); print(c(za, zb))
print(RTMB::pnorm(c(za, zb), log.p = TRUE)); print(RTMB::pnorm(c(za,zb), lower.tail=FALSE, log.p=TRUE))
print(RTMB::logspace_add(-3, -Inf)); print(E$log1m_exp_ad(c(-1, 0, 2)))
