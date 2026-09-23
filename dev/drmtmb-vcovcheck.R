# vcov(fit, full = TRUE) is frmtmb's covariance of the outer parameter
# vector in opt$par order: checked against RTMB::sdreport() directly,
# because the agreement engine reads standard errors from it.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
d <- sim_grouped(101)
chk <- function(ff) {
  a <- sqrt(diag(vcov(ff, full = TRUE)))
  s <- RTMB::sdreport(ff$obj, par.fixed = ff$opt$par)
  b <- sqrt(diag(s$cov.fixed))
  cat("outer parameters:", paste(names(ff$opt$par), collapse = " "), "\n")
  cat("max |log(se_vcov / se_sdreport)|:",
      format(max(abs(log(unname(a) / unname(b)))), digits = 3), "\n")
}
chk(frm(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d))
chk(frm(bf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1), family = student(),
        data = d))
