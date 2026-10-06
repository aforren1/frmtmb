# Punch round 1: scope of the empty-fixed-part fitted() failure. Is it
# gp() only, or any predictor whose only terms are penalized/grouped?
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(11)
n <- 120
d <- data.frame(x = round(runif(n, 0, 5), 1), z = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)))
lat <- sin(d$x) + 0.5 * d$z + rlogis(n)
d$y <- factor(cut(lat, quantile(lat, c(0, 0.3, 0.6, 1)),
                  include.lowest = TRUE, labels = FALSE), ordered = TRUE)
d$yg <- sin(d$x) + rnorm(n, 0, 0.3)
try1 <- function(lab, f, fam) {
  fit <- tryCatch(suppressWarnings(frm(f, data = d, family = fam)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(lab, "| fit error:", conditionMessage(fit), "\n")
    return(invisible())
  }
  r <- tryCatch(fitted(fit), error = function(e) e)
  cat(lab, "| fitted():", if (inherits(r, "error")) conditionMessage(r)
      else "ok", "\n")
}
try1("cumulative disc ~ 0 + gp(x, k = 6)",
     bf(y ~ z, disc ~ 0 + gp(x, k = 6)), cumulative())
try1("cumulative disc ~ 0 + s(x, k = 5)",
     bf(y ~ z, disc ~ 0 + s(x, k = 5)), cumulative())
try1("cumulative disc ~ 0 + (1 | g)",
     bf(y ~ z, disc ~ 0 + (1 | g)), cumulative())
try1("gaussian sigma ~ 0 + gp(x, k = 6)",
     bf(yg ~ z, sigma ~ 0 + gp(x, k = 6)), gaussian())
try1("gaussian sigma ~ 0 + (1 | g)",
     bf(yg ~ z, sigma ~ 0 + (1 | g)), gaussian())
try1("gaussian yg ~ 0 + gp(x, k = 6)", bf(yg ~ 0 + gp(x, k = 6)),
     gaussian())
