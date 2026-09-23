# Lane wt-reunc: fitted() on a mixed ordinal fit. The category
# probabilities take the finite-difference delta method; this checks
# the new b columns against a brute-force Monte Carlo over the SAME
# joint law, and times the route.
#
#   Rscript dev/reunc-ordinal.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(5)
G <- 10; m <- 12
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
lat <- 0.8 * d$x + rnorm(G, 0, 1)[d$g] + rlogis(G * m)
d$y <- cut(lat, c(-Inf, -1, 0.5, Inf), labels = FALSE)
d$y <- factor(d$y, ordered = TRUE)
fit <- frm(bf(y ~ x + (1 | g)) + cumulative(), data = d)
nd <- data.frame(x = c(0, 1), g = factor(c(1, 2), levels = levels(d$g)))
t0 <- proc.time()[["elapsed"]]
f <- fitted(fit, newdata = nd)
t1 <- proc.time()[["elapsed"]]
fna <- fitted(fit, newdata = nd, re_formula = NA)
cat("fitted() seconds:", t1 - t0, "\n")
cat("Est.Error at known level:\n"); print(f[, "Est.Error", ])
cat("Est.Error re_formula = NA:\n"); print(fna[, "Est.Error", ])
# brute force: draw (outer, b) from the joint covariance the delta
# method uses and push each draw through fitted()'s own estimate
ds <- frmtmb:::fit_draw_space(fit)
map <- ds$map
jc <- frmtmb:::get_joint_cov(fit)
rn <- jc$names
pos <- c(unlist(lapply(unique(map$comp), function(cp) which(rn == cp))),
         which(rn == "b"))
V <- as.matrix(jc$V[pos, pos])
v0 <- c(frmtmb:::fit_outer_vector(fit, map), fit$estimates[["b"]])
p <- length(map$names)
L <- chol(V)
N <- 4000
set.seed(7)
draws <- matrix(NA_real_, N, length(f[, 1, ]))
for (s in seq_len(N)) {
  v <- v0 + drop(crossprod(L, rnorm(length(v0))))
  fs <- frmtmb:::fit_set_outer(fit, v[seq_len(p)], map)
  fs$estimates[["b"]][] <- v[-seq_len(p)]
  draws[s, ] <- as.vector(frmtmb:::fitted_point(fs, nd))
}
mc <- matrix(apply(draws, 2, sd), nrow(nd))
cat("Monte Carlo sd over the joint law:\n"); print(round(mc, 5))
cat("ratio delta / MC:\n"); print(round(f[, "Est.Error", ] / mc, 4))
