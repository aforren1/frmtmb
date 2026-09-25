# Localize the NaN gradient the first hunt found at a = exp(10).
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib2",
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE,
                                   export_all = TRUE))
q <- c(0.3, 0.45, 0.7, 1.0)
tt <- function(p) ddm_floor(q - p[3], ddm_gng_t_floor)
fs <- list(
  small = function(p) sum(ddm_cdf_small(tt(p), p[1], exp(p[2]), 0.5)),
  small_lo = function(p) sum(ddm_lower_cdf_small(tt(p), p[1], exp(p[2]), 0.5)),
  small_up = function(p) sum(ddm_lower_cdf_small(tt(p), -p[1], exp(p[2]), 0.5)),
  large = function(p) sum(ddm_surv_large(tt(p), p[1], exp(p[2]), 0.5)),
  logu = function(p) sum(log(ddm_floor(tt(p) / exp(2 * p[2]), ddm_u_floor))))
for (p in list(c(-60, 10, 0.3), c(-60, 10, 0.29999999), c(-60, 10, 0.5), c(-5, 10, 0.3), c(-60, 3, 0.3))) {
  for (nm in names(fs)) {
    tp <- RTMB::MakeTape(fs[[nm]], c(0, 0, 0))
    cat(format(p), nm, format(tp(p), digits = 4), "|",
        format(tp$jacobian(p), digits = 4), "\n")
  }
}
