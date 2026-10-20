# Lane illse2, item 2: why the route does not take degree 7 seed 4
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
set.seed(4); n <- 300; deg <- 7
d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
xn <- paste0("x", seq_len(deg))
for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
f <- frm(fo, data = d)
p <- f$opt$par
h <- ns$fit_outer_hessian(f)
H <- h$H
cat("se_explained:", format(f$cache$se_explained), "| analysis:",
    !is.null(f$cache$se_analysis), "| hessian_fixed clean:",
    format(f$cache$hessian_fixed$clean), "\n")
cat("edge:", ns$se_edge_sd(f, H), "| held:", ns$se_bound_held(f, p), "\n")
cat("could_act:", ns$se_tier3_could_act(f, H, h$E, FALSE,
    ns$se_bound_held(f, p), ns$se_edge_sd(f, H)), "\n")
M <- ns$reparam_basis(f$frame, names(p))
cat("defer:", ns$se_reparam_defer(f, H, attr(M, "pos")), "\n")
D <- sqrt(abs(diag(H)))
ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
cat("unit-diag eigen range:", signif(range(ev), 3), "\n")
rp <- ns$se_reparam_cov(f, H, p, defer = FALSE)
cat("route without defer:", !is.null(rp), "\n")
an <- ns$cov_from_hessian(f, H, h$E)
cat("cov_from_hessian tier:", an$tier, "root:", !is.null(an$root), "\n")
