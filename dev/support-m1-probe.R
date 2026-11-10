# Lane support, review finding 2: evaluations on doubles that feed a fit
# must keep the fitted objective's arithmetic. Two direct cases:
#  1. check_custom_family() on the GEV at test values with a row outside
#     the support (R/interop.R:74): its plain-numeric value must be the
#     taped one.
#  2. the Frechet log density where nu z < -700 (the exp_neg_capped()
#     cap), on doubles against the tape: the value nl_domain_start()
#     and nl_coef_curved() read through nll() on doubles.
# Then the same on doubles with a reader's exact = TRUE, which must
# differ there.
#
#   Rscript dev/support-m1-probe.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
gev <- gen_extreme_value()
y <- c(0, 0.5, 2)
r <- tryCatch({
  check_custom_family(gev, y = y,
                      dpars = list(mu = c(1, 1, 1), sigma = c(0.1, 1, 1),
                                   xi = c(0.3, 0.3, 0.3)))
  "passes"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("check_custom_family, GEV row 1 outside:", r, "\n")
fr <- frechet()
# w = nu z with z = log(y / scale): y tiny puts w near -705
dp <- list(mu = 1, nu = 3)
s <- 1 / gamma(1 - 1 / 3)
yy <- s * exp(-705 / 3)
d1 <- fr$lpdf(yy, dp, list())
tp <- RTMB::MakeTape(function(p) fr$lpdf(yy, list(mu = p[1], nu = p[2]),
                                         list()), c(1, 3))
cat(sprintf("frechet lpdf at w = -705: doubles %.6g, tape %.6g\n", d1,
            tp(c(1, 3))))
has_exact <- "exact" %in% names(formals(frmtmb::row_lpdf))
if (has_exact) {
  ex <- frmtmb::row_lpdf(fr, yy, yy, dp, list(), NULL, exact = TRUE)
  cat(sprintf("frechet reader (exact = TRUE): %.6g\n", ex))
}
