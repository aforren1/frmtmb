# Reviewer, fix round 1: is the compois objective rough at nu = e^3 and
# e^4 because of the censored rows, or already in the density? Central
# differences at three steps, the data of cdfs-rev-grad.R (group 3).
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-cdfs/dev/cdfs-rev-grad.R")
Sys.setenv(GROUP = 3)
eval(parse(text = src[1:(grep("^pick <- ", src) - 1)]))
d <- fams$compois$d
for (keep in c("none", "all", "-1", "1", "2")) {
  dd <- d
  if (keep == "none") dd$cc <- 0 else if (keep != "all") dd$cc[dd$cc != as.numeric(keep)] <- 0
  ob <- frm(bf(y | cens(cc, y2) ~ x, nu ~ 1), family = brmsfamily("com_poisson"), data = dd,
            dry_run = "objective")$obj
  for (v in c(3, 4)) {
    p <- ob$par; p[3] <- v
    g <- as.numeric(ob$gr(p))[3]
    fds <- vapply(c(1e-3, 1e-5), function(h) { e <- c(0, 0, h); (ob$fn(p + e) - ob$fn(p - e)) / (2 * h) }, 0)
    cat(sprintf("codes %-4s nu eta %g: fn %.8g  d/d eta_nu AD %.6g FD(1e-3) %.6g FD(1e-5) %.6g\n",
                keep, v, ob$fn(p), g, fds[1], fds[2]))
  }
}
ns <- asNamespace("frmtmb")
cat("== RTMB::dcompois2 and compois_cdf_log at nu = 20 and 55, mu 2.7, along nu\n")
for (nu in c(20, 20 * (1 + 1e-5), 55, 55 * (1 + 1e-5))) {
  cat(sprintf("nu %.6f: lpmf(0:4) %s | lF(2) %.12g lS(2) %.12g\n", nu,
              paste(signif(RTMB::dcompois2(0:4, 2.7, nu, log = TRUE), 10), collapse = " "),
              ns$compois_cdf_log(2, 2.7, nu)$lF, ns$compois_cdf_log(2, 2.7, nu)$lS))
}
