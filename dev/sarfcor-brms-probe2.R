# brms 2.23.0 and new data for sar() and fcor(): an unfitted brmsfit
# (empty = TRUE, no compile) read by standata(newdata =, newdata2 =).
# Run: Rscript dev/sarfcor-brms-probe2.R > dev/sarfcor-out/brms-probe2.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 9
xy <- expand.grid(i = 1:3, j = 1:3)
A <- (as.matrix(dist(xy)) == 1) * 1
W <- A / rowSums(A)
V <- 0.5^abs(outer(1:n, 1:n, "-"))
d <- data.frame(y = rnorm(n), x = rnorm(n))
try_ <- function(lab, expr) {
  out <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("---", lab, "\n")
  if (is.character(out)) cat(out, sep = "\n") else print(out)
}
fs <- brm(y ~ x + sar(W), data = d, data2 = list(W = W), empty = TRUE)
ff <- brm(y ~ x + fcor(V), data = d, data2 = list(V = V), empty = TRUE)
nd <- d[1:4, ]
try_("sar newdata, no newdata2", dim(standata(fs, newdata = nd)$Msar))
try_("sar newdata, newdata2 W[1:4, 1:4]",
     dim(standata(fs, newdata = nd,
                  newdata2 = list(W = W[1:4, 1:4]))$Msar))
try_("fcor newdata, no newdata2", dim(standata(ff, newdata = nd)$Mfcor))
try_("fcor newdata, newdata2 V[1:4, 1:4]",
     dim(standata(ff, newdata = nd,
                  newdata2 = list(V = V[1:4, 1:4]))$Mfcor))
# the NA rows: which rows of M brms keeps
dn <- d
dn$x[c(2, 5)] <- NA
st <- suppressWarnings(standata(y ~ x + fcor(V), data = dn,
                                data2 = list(V = V)))
cat("--- fcor with x[c(2, 5)] NA: Mfcor equals V[-c(2, 5), -c(2, 5)]:",
    isTRUE(all.equal(unname(st$Mfcor), V[-c(2, 5), -c(2, 5)])), "\n")
st <- suppressWarnings(standata(y ~ x + sar(W), data = dn,
                                data2 = list(W = W)))
cat("--- sar with x[c(2, 5)] NA: Msar equals W[-c(2, 5), -c(2, 5)]:",
    isTRUE(all.equal(unname(st$Msar), unname(W[-c(2, 5), -c(2, 5)]))),
    "; row sums after the drop:", format(rowSums(st$Msar), digits = 3),
    "\n")
# variables brms would name: the parameter block of the Stan code
for (f in list(y ~ x + sar(W), y ~ x + sar(W, type = "error"),
               y ~ x + fcor(V))) {
  code <- stancode(f, data = d, data2 = list(W = W, V = V))
  pb <- regmatches(code, regexpr("parameters \\{[^}]*\\}", code))
  cat("--- parameters of", deparse(f), "\n", pb[1], "\n")
}
# sigma with a covariance M: brms's Cov = sigma^2 * M
code <- stancode(y ~ x + fcor(V), data = d, data2 = list(V = 4 * V))
cat("--- fcor likelihood line:",
    grep("fcor_hom_lpdf\\(Y", strsplit(code, "\n")[[1]], value = TRUE), "\n")
