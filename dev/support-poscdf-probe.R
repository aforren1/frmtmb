# Lane support, survey row: a trunc() lower bound below 0 on a family with
# a positive response. P(Y > lb) is 1 there, so trunc(lb = -1) is no
# truncation at all and the fit must equal the untruncated one. Seed 401,
# n 300, y ~ exponential with mean 2.
#
#   Rscript dev/support-poscdf-probe.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
set.seed(401)
n <- 300
d <- data.frame(y = rexp(n, 1 / 2), lb = -1, cy = 0.5)
fams <- list(exponential = exponential(), gamma = Gamma(link = "log"),
             weibull = weibull(), lognormal = lognormal(),
             frechet = frechet(), inverse.gaussian = inverse.gaussian())
for (nm in names(fams)) {
  f0 <- tryCatch(frm(y ~ 1, family = fams[[nm]], data = d),
                 error = function(e) conditionMessage(e))
  f1 <- tryCatch(suppressWarnings(
    frm(y | trunc(lb = lb) ~ 1, family = fams[[nm]], data = d)),
    error = function(e) conditionMessage(e))
  if (is.character(f0) || is.character(f1)) {
    cat(sprintf("%-16s untruncated: %s | trunc(lb = -1): %s\n", nm,
                if (is.character(f0)) f0 else "fits",
                if (is.character(f1)) substr(f1, 1, 120) else "fits"))
    next
  }
  cat(sprintf(paste("%-16s logLik untruncated %.10g, trunc(lb = -1) %.10g;",
                    "Intercept %.8g vs %.8g\n"), nm,
              as.numeric(logLik(f0)), as.numeric(logLik(f1)),
              fixef(f0)[1, 1], fixef(f1)[1, 1]))
}
# the CDF slots themselves at q = -1 and q = 0, on doubles
for (nm in c("exponential", "weibull", "frechet")) {
  fam <- fams[[nm]]
  dp <- switch(nm, exponential = list(mu = 2),
               weibull = list(mu = 2, shape = 1.5),
               frechet = list(mu = 2, nu = 3))
  F <- tryCatch(fam$lcdf(c(-1, 0), dp, list()),
                error = function(e) conditionMessage(e))
  lS <- tryCatch(fam$lccdf(c(-1, 0), dp, list()),
                 error = function(e) conditionMessage(e))
  cat(sprintf("%-12s F(-1), F(0): %s; log S(-1), log S(0): %s\n", nm,
              paste(format(F, digits = 6), collapse = ", "),
              paste(format(lS, digits = 6), collapse = ", ")))
}
