# Lane wt-resmooth. The PRE-EXISTING half of the finite-difference
# standard-error gap: on an ordinal fit with a POPULATION smooth, the
# category-probability Est.Error differenced no smooth coefficients at
# ANY re_formula, so it left out the curve's uncertainty while the
# analytic route carried it. smooth_b_idx() closes this too.
#   RESMOOTH_LIB=base Rscript dev/resmooth-popse.R > dev/resmooth-popse-before.txt
#   Rscript dev/resmooth-popse.R > dev/resmooth-popse-after.txt
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")
set.seed(23)
n <- 240
d <- data.frame(x = stats::runif(n, -2, 2))
lat <- 1.2 * sin(2 * d$x) + stats::rlogis(n)
d$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
              ordered = TRUE)
fit <- suppressWarnings(frm(bf(y ~ s(x, k = 8)) + cumulative(), data = d))
nd <- data.frame(x = c(-1.5, 0, 1.5))
for (rf in list(NULL, NA)) {
  lab <- if (is.null(rf)) "NULL" else "NA"
  got <- fitted(fit, newdata = nd, re_formula = rf)[, "Est.Error", ]
  f <- function(x) frmtmb:::fitted_point(x, nd, re_formula = rf)
  bl <- fit$frame[["re_blocks"]]
  sm <- sort(unique(unlist(lapply(bl[vapply(bl, function(b) {
    b[["covstruct"]] %in% c("smooth", "gp", "hsgp")
  }, NA)], `[[`, "b_idx"))))
  ref <- as.vector(frmtmb:::fit_fd_se(fit, f, b_idx = sm, b_batch = NULL))
  rel <- (ref - as.vector(got)) / ref
  cat(sprintf(paste0("re_formula = %-4s | shipped Est.Error %s | ",
                     "with the smooth's %d coefficients %s | ",
                     "median shortfall %.4f\n"),
              lab, paste(sprintf("%.5f", as.vector(got)), collapse = " "),
              length(sm),
              paste(sprintf("%.5f", ref), collapse = " "),
              stats::median(rel)))
}
