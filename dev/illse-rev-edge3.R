# Reviewer of lane illse: deg 8 seed 3 of dev/illse-rev-edge2.R, the sd at
# 8e-23 that took the route: the SE of log sd and what the edge test saw.
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
deg <- 8; seed <- 3
set.seed(seed)
n <- 300
d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
xn <- paste0("x", seq_len(deg))
for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
r <- cap(frm(fo, data = d))
f <- r$v
print(f$opt$par)
print(f$opt$polish)
V <- suppressWarnings(vcov(f, full = TRUE))
print(sqrt(diag(V)))
print(r$w)
s <- cap(summary(f)); print(s$v$random); print(s$w)
vc <- cap(VarCorr(f)); print(vc$v); print(vc$w)
g <- f$obj$gr(f$opt$par); print(g)
H <- stats::optimHess(f$opt$par, f$obj$fn, f$obj$gr)
print(diag(H))
if (exists("se_edge_sd", ns)) print(ns$se_edge_sd(f, H))
