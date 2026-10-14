# Lane illse, prototype debug: M' H M by differences against the exact
# one (256-bit) on degree 8, seed 1, and the gradient against the exact
# gradient at the fit's point.
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
ns <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-illse/R/se-reparam.R",
           envir = ns)
deg <- as.integer(Sys.getenv("DEG", "8"))
ds <- des_poly(deg)(1)
f <- suppressWarnings(frm(ds$formula, data = ds$data))
p <- f$opt$par
X <- as.matrix(f$frame$linpreds[[1]]$X)
y <- ds$data$y
b <- p[names(p) == "beta"]; ls <- p[names(p) == "betad"]
Xm <- mpfrArray(X, mp_prec, dim = dim(X))
r <- mp(y) - Xm %*% mp(b)
s2 <- exp(2 * mp(ls))
ge <- c(asNumeric(-crossprod(Xm, r) / s2), asNumeric(length(y) - sum(r^2) / s2))
g <- f$obj$gr(p)
cat("gradient frmtmb vs exact:\n"); print(rbind(frmtmb = drop(g), exact = ge))
cat("objective frmtmb", f$opt$objective, " exact nll ",
    asNumeric(length(y) * mp(ls) + sum(r^2) / (2 * s2) +
              length(y) / 2 * log(2 * pi)), "\n")
M <- ns$reparam_basis(f$frame, names(p))
h <- ns$reparam_hessian(f$obj, p, M)
np <- length(p)
He <- mpfrArray(0, mp_prec, dim = c(np, np))
He[1:(np - 1), 1:(np - 1)] <- crossprod(Xm) / s2
hb <- -2 * crossprod(Xm, r) / s2
He[1:(np - 1), np] <- hb
He[np, 1:(np - 1)] <- hb
He[np, np] <- 2 * sum(r^2) / s2
Mm <- mpfrArray(M, mp_prec, dim = dim(M))
Hte <- asNumeric(t(Mm) %*% He %*% Mm)
cat("M'HM by differences, last column (log sigma):\n")
print(rbind(fd = h$H[, np], exact = Hte[, np]))
cat("diag:\n"); print(rbind(fd = diag(h$H), exact = diag(Hte)))
