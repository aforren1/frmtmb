# Review: the `theta thetaac_1` row of default_prior() for a sar() and
# an ar() model, on the library given as argument 1 (then the user lib).
a <- commandArgs(TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
has_sp <- suppressWarnings(requireNamespace("frmtmb.spatial", quietly = TRUE))
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "spatial:", has_sp, "\n")
ptab <- function(p) as.data.frame(p)[, c("prior", "class", "coef", "lb",
                                         "ub", "source")]
k <- 5
g <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A)
set.seed(1)
d <- data.frame(x = rnorm(25), y = rnorm(25), t = 1:25, s = 1)
if (has_sp) {
  suppressMessages(library(frmtmb.spatial))
  print(ptab(default_prior(y ~ x + sar(W), data = d, data2 = list(W = W))))
  fit <- frm(y ~ x + sar(W), data = d, data2 = list(W = W))
  cat("prior on theta coef thetaac_1:\n")
  r <- tryCatch(validate_prior(set_prior("normal(0, 1)", class = "theta",
                                         coef = "thetaac_1"),
                               y ~ x + sar(W), data = d, data2 = list(W = W)),
                error = function(e) conditionMessage(e))
  print(r)
}
cat("ar():\n")
print(ptab(default_prior(y ~ x + ar(t, s), data = d)))
V <- 0.5^abs(outer(1:25, 1:25, "-"))
if ("fcor" %in% getNamespaceExports("frmtmb")) {
  cat("fcor():\n")
  print(ptab(default_prior(y ~ x + fcor(V), data = d, data2 = list(V = V))))
}
