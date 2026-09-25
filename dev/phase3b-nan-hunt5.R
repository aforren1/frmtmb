# Localize the non-finite values phase3b-nan-hunt4.R found.
a0 <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a0[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.eam")
chk <- function(label, f, p) {
  tp <- RTMB::MakeTape(f, p)
  cat(sprintf("%-28s value %-12s grad %s\n", label, format(tp(p), digits = 5),
              paste(format(tp$jacobian(p), digits = 4), collapse = " ")))
}
pts <- list(c(-60, 6, 0.45), c(-60, 6, -5), c(60, 6, 0.45), c(-20, 3, 0.45),
            c(-5, 6, -5), c(0, 6, -5))
for (p in pts) {
  cat("== v", p[1], "log a", p[2], "t0", p[3], "\n")
  tt <- function(p) ns$ddm_floor(0.45 - p[3], ns$ddm_gng_t_floor)
  chk("lsurv_small", function(p) ns$ddm_lsurv_small(tt(p), p[1], exp(p[2]), 0.5), p)
  chk("lsurv_large", function(p) ns$ddm_lsurv_large(tt(p), p[1], exp(p[2]), 0.5), p)
  chk("llower_small", function(p) ns$ddm_llower_small(tt(p), p[1], exp(p[2]), 0.5), p)
  chk("llower_large", function(p) ns$ddm_llower_large(tt(p), p[1], exp(p[2]), 0.5), p)
  chk("llower_prob", function(p) ns$ddm_llower_prob(p[1], exp(p[2]), 0.5), p)
}
