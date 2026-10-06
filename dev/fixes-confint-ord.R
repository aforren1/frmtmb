# Lane fixes, item 4c: what confint() reports for an ordinal fit's
# thresholds, beside fixef() and vcov(full = TRUE).
#   Rscript dev/fixes-confint-ord.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
fits <- list(
  cum = frm(y ~ x, family = cumulative(), data = d),
  sratio = frm(y ~ x, family = sratio(), data = d),
  cum_equi = frm(y ~ x, family = cumulative(threshold = "equidistant"),
                 data = d),
  acat_equi = frm(y ~ x, family = acat(threshold = "equidistant"), data = d),
  cum_stz = frm(y ~ x, family = cumulative(threshold = "sum_to_zero"),
                data = d),
  cum_gr = frm(yg | thres(gr = g) ~ x, family = cumulative(), data = d),
  cum_disc = frm(bf(y ~ x, disc ~ z), family = cumulative(), data = d))
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("==", nm, "\n")
  print(round(confint(f), 5))
  cat("fixef:\n"); print(round(fixef(f), 5))
  cat("vcov(full = TRUE) names:", rownames(vcov(f, full = TRUE)), "\n")
  for (p in c("Intercept[1]", "Intercept[2]", "b_Intercept[1]", "delta")) {
    r <- tryCatch(rownames(suppressMessages(confint(f, parm = p))),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("  confint(parm = %s): %s\n", p, substr(r, 1, 150)))
  }
}
