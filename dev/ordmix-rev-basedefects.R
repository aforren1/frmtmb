# Reviewer of lane ordmix: the two base defects the lane fixed, on the
# arm given, and every other route to the cs() check it changed.
# Usage: Rscript dev/ordmix-rev-basedefects.R <lane|base>. Seed 20261097.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261097)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
lat <- 0.8 * d$x + rlogis(n)
d$y <- 1L + (lat > -1) + (lat > 0) + (lat > 1)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
d$yg <- d$w + rnorm(n)
try_fit <- function(lab, expr) {
  r <- tryCatch({
    f <- suppressWarnings(expr)
    sprintf("FITS logLik %.10f; coefs %s", as.numeric(logLik(f)),
            paste(rownames(fixef(f)), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", lab, substr(r, 1, 260)))
}
try_fit("sratio disc ~ cs(z)",
        frm(bf(y ~ x, disc ~ cs(z)), family = sratio(), data = d))
try_fit("sratio y ~ x + cs(z)",
        frm(bf(y ~ x + cs(z)), family = sratio(), data = d))
try_fit("sratio y ~ x",
        frm(bf(y ~ x), family = sratio(), data = d))
try_fit("gaussian sigma ~ cs(w)",
        frm(bf(yg ~ w, sigma ~ cs(w)), family = gaussian(), data = d))
try_fit("gaussian mu ~ cs(w)",
        frm(bf(yg ~ cs(w)), family = gaussian(), data = d))
try_fit("sratio nlpar a ~ cs(z)",
        frm(bf(y ~ a, a ~ 0 + cs(z), nl = TRUE), family = sratio(),
            data = d))
try_fit("hurdle hu ~ cs(z)",
        frm(bf(yh ~ x, hu ~ cs(z)), family = hurdle_cumulative(), data = d))
try_fit("acat disc ~ 0 + cs(z)",
        frm(bf(y ~ x, disc ~ 0 + cs(z)), family = acat(), data = d))
try_fit("cumulative y ~ cs(x)",
        frm(bf(y ~ cs(x)), family = cumulative(), data = d))
try_fit("cratio y ~ cs(x) (control)",
        frm(bf(y ~ cs(x)), family = cratio(), data = d))

## defect 2: brms_par_labels() with no location column
f2 <- suppressWarnings(frm(bf(y ~ 1, disc ~ 0 + z), family = cumulative(),
                           data = d))
lab <- frmtmb:::brms_par_labels(f2)
cat("labels y ~ 1, disc ~ 0 + z:", paste(lab, collapse = " "), "\n")
cat("  opt$par names:", paste(names(f2$opt$par), collapse = " "), "\n")
f3 <- suppressWarnings(frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(),
                           data = d))
cat("labels y ~ x, disc ~ 0 + z:",
    paste(frmtmb:::brms_par_labels(f3), collapse = " "), "\n")
