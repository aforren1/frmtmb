# Reviewer reconnaissance: what does fitted() return on an ordinal fit,
# and how long does one fitted_point() call cost (for sizing the Monte
# Carlo reference)?
LIB <- if (identical(Sys.getenv("REVLIB"), "base")) {
  "C:/Users/adf44/source/r/rellib-r3"
} else "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")

set.seed(23)
n <- 240
d1 <- data.frame(x = stats::runif(n, -2, 2))
lat <- 1.2 * sin(2 * d1$x) + stats::rlogis(n)
d1$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
f1 <- suppressWarnings(frm(bf(y ~ s(x, k = 8)), family = cumulative(),
                          data = d1))
nd <- data.frame(x = c(-1.5, 0, 1.5))
a <- fitted(f1, newdata = nd, re_formula = NA)
cat("class:", class(a), "dim:", paste(dim(a), collapse = "x"), "\n")
print(dimnames(a))
print(a)
cat("--- args of fitted method ---\n")
print(names(formals(frmtmb:::fitted.frmtmb_fit)))
t0 <- proc.time()
for (i in 1:20) invisible(frmtmb:::fitted_point(f1, nd, re_formula = NA))
cat("20 fitted_point calls, seconds:", (proc.time() - t0)[["elapsed"]], "\n")
t0 <- proc.time()
s1 <- fitted(f1, newdata = nd, re_formula = NA)
cat("one fitted() with Est.Error, seconds:",
    (proc.time() - t0)[["elapsed"]], "\n")
cat("smooth_b_idx exists:", exists("smooth_b_idx", asNamespace("frmtmb")),
    "\n")
