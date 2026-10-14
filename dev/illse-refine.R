# Lane illse: the first pass of reparam_hessian_refined() one Newton
# step from the stop of a cumulative degree-8 fit (seed 1): the
# threshold direction is below its own noise there.
#   Rscript dev/illse-refine.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
set.seed(1)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- factor(cut(2 * sin(3 * d$x) + stats::rlogis(200),
                  c(-Inf, -1, 0, 1, Inf), labels = FALSE), ordered = TRUE)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
f <- suppressWarnings(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = d,
             family = cumulative()))
ns <- asNamespace("frmtmb")
p <- f$opt$par; obj <- f$obj
M <- ns$reparam_basis(f$frame, names(p))
h <- ns$reparam_hessian_refined(obj, p, M)
g <- drop(obj$gr(p)); gs <- drop(crossprod(h$M, g))
R <- chol(h$H); st <- -drop(backsolve(R, forwardsolve(t(R), gs)))
q <- p + drop(h$M %*% st)
h1 <- ns$reparam_hessian(obj, q, M)
D <- sqrt(diag(h1$H)); e <- eigen(h1$H / outer(D, D), symmetric = TRUE)
cat("first pass at the stepped point, eigenvalues / largest:\n")
print(e$values / max(e$values))
nk <- sqrt(colSums(((h1$E / outer(D, D)) %*% e$vectors)^2))
cat("their noise:\n"); print(nk)
