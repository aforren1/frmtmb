# Lane seflags, item 5: ls ~ mo(inc) * age, seed 2 of dev/optima-segap.R,
# whose interaction simplex sits at a face: the three smallest
# directions, their gradient curvature at steps 1e-3, 1e-2, 1e-1 and
# se_ridge_curvature()'s corrected value, and the simplex summary.
#   Rscript dev/seflags-segapface.R
.libPaths(c("C:/Users/adf44/source/r/wt-seflags-lib", "C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); ns <- asNamespace("frmtmb")
set.seed(2); inc <- sample(c(0L, 1L, 3L), 100, TRUE)
d <- data.frame(inc, ls = c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7), age = rnorm(100, 40, 10))
f <- suppressWarnings(frm(ls ~ mo(inc) * age, data = d))
p <- f$opt$par; print(round(p, 4)); cat("grad", signif(f$obj$gr(p), 2), "\n")
H <- f$obj$he(p); H <- (H + t(H)) / 2
D <- sqrt(abs(diag(H))); e <- eigen(H / outer(D, D), symmetric = TRUE)
big <- max(abs(e$values)); cat("ev rel", signif(e$values / big, 3), "\n")
np <- length(p)
for (k in (np - 2):np) {
  vk <- e$vectors[, k]
  cur <- function(z, h) { dd <- h * vk / D
    sum((as.numeric(f$obj$gr(z + dd)) - as.numeric(f$obj$gr(z - dd))) / D * vk) / (2 * h) }
  cat("k", k, "ev", signif(e$values[k], 3), "v", signif(vk, 2), "\n  curv:", signif(sapply(c(1e-3, 1e-2, 1e-1), function(h) cur(p, h)), 3),
      "c0", signif(ns$se_ridge_curvature(f, p, seq_len(np), e, D, k, rep(TRUE, np)), 3), "\n")
}
print(suppressWarnings(summary(f))$mo)
