# Lane illse2, item 4: why a bound-held raw polynomial stops short
# (dev/illse2-p4bound.R, degree 8 seed 1 and degree 7 seed 3).
#   ILLSE2_ARM=lane Rscript dev/illse2-p4bdbg.R <deg> <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
deg <- as.integer(a[1]); seed <- as.integer(a[2])
set.seed(seed)
d <- data.frame(x = runif(200, 1, 2))
d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
xn <- paste0("x", 1:deg)
for (k in 1:deg) d[[xn[k]]] <- d$x^k
top <- xn[deg]
fo <- as.formula(paste("y ~", paste0(xn, collapse = " + ")))
l0 <- lm(fo, data = d)
ub <- coef(l0)[[top]] - sqrt(diag(vcov(l0)))[[top]]
pr <- set_prior("", class = "b", coef = top, ub = ub)
f <- frm(fo, data = d, prior = pr,
         control = frmtmb_control(verbose = TRUE))
p <- f$opt$par
cat("polish:", format(f$opt$polish), "\n")
g <- f$obj$gr(p)
cat("gradient:", signif(g, 3), "\n")
box <- ns$fit_outer_box(f)
cat("active:", ns$grad_bound_active(p, g, box), "\n")
H <- f$obj$he(p)
cat("hess_needs_reparam (all):", ns$hess_needs_reparam(H), "\n")
fr <- seq_along(p)[-(deg + 1)]
cat("hess_needs_reparam (free):", ns$hess_needs_reparam(H[fr, fr]), "\n")
D <- sqrt(abs(diag(H)))
ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
cat("unit-diag eigen range:", signif(range(ev), 3), "\n")
H2 <- (H + t(H)) / 2
hh <- ns$headroom_reparam_hessian(f, H2, seq_along(p))
cat("headroom hessian: H", dim(hh$H), "M", dim(hh$M), "
")
print(traceback)
cat("conv_headroom:", tryCatch(ns$conv_headroom(f, g), error = function(e) conditionMessage(e)), "| grad_headroom:",
    ns$grad_headroom(f, g, ns$grad_bound_active(p, g, box) %in% TRUE),
    "\n")
