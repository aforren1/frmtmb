# Lane seflags: dev/optima-segap.R's ls ~ mo(inc) + age, seed 3, which
# stays silent: its spectrum, the corrected curvature of the smallest
# direction, the gradient and what the exact two-step move costs.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); ns <- asNamespace("frmtmb")
s <- if (length(a) > 1) as.integer(a[2]) else 3L
set.seed(s)
inc <- sample(c(0L, 1L, 3L), 100, TRUE)
ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
d <- data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
f <- suppressWarnings(frm(ls ~ mo(inc) + age, data = d))
p <- f$opt$par; print(p); cat("grad", signif(f$obj$gr(p), 3), "\n")
H <- f$obj$he(p); H <- (H + t(H)) / 2
D <- sqrt(abs(diag(H))); e <- eigen(H / outer(D, D), symmetric = TRUE)
big <- max(abs(e$values)); cat("ev rel", signif(e$values / big, 3), "\n")
k <- length(p)
cat("v:", signif(e$vectors[, k], 3), "\n")
c0 <- ns$se_ridge_curvature(f, p, seq_len(k), e, D, k, rep(TRUE, k))
cat("c0", c0, "thr", 1e-9 * big, "\n")
print(ns$mo_simplex_weights %||% NULL)
print(suppressWarnings(summary(f))$mo)
vk <- e$vectors[, k]
cur <- function(z, h) { dd <- h * vk / D
  sum((as.numeric(f$obj$gr(z + dd)) - as.numeric(f$obj$gr(z - dd))) / D * vk) / (2 * h) }
cat("curv at p:", signif(sapply(c(1e-4, 1e-3, 1e-2, 1e-1, 1), function(h) cur(p, h)), 3), "\n")
f0 <- f$obj$fn(p)
cat("f along v (unit diag steps 1e-2..10):", signif(sapply(c(1e-2, 0.1, 1, 10), function(h) f$obj$fn(p + h * vk / D) + f$obj$fn(p - h * vk / D) - 2 * f0), 3), "\n")
cat("zeta D:", signif(D[5:6], 3), "\n")
