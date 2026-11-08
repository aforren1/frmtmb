# Lane bigx, step 2: the intermediates of one poisson fit at x + m.
#   BIGX_ARM=before Rscript dev/bigx-trace.R [m] [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
a <- commandArgs(TRUE)
m <- if (length(a)) as.numeric(a[1]) else 1e6
s <- if (length(a) > 1) as.integer(a[2]) else 2L
d <- bigx_data("poisson", m, s)
fc <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
pc <- fc$opt$par
# the centered optimum in the raw coordinates
pr <- pc
pr[1] <- pc[1] - pc[2] * m
cat("centered optimum:", format(pc, digits = 10), "\n")
cat("as raw:          ", format(pr, digits = 10), "\n")
f <- suppressWarnings(frm(bf(y ~ w + (1 | g)), family = poisson(),
                          data = d, control = frmtmb_control(verbose = 2)))
cat("raw stop:        ", format(f$opt$par, digits = 10), "\n")
cat("objective: raw stop", format(f$opt$objective, digits = 12),
    " raw tape at centered optimum", format(f$obj$fn(pr), digits = 12),
    " centered", format(fc$opt$objective, digits = 12), "\n")
cat("start (frame template beta):", format(f$start %||% NA), "\n")
cat("par_units:", format(f$par_units %||% 1), "\n")
cat("gradient at stop:", format(as.numeric(f$obj$gr(f$opt$par))), "\n")
cat("gradient at centered optimum (raw tape):",
    format(as.numeric(f$obj$gr(pr))), "\n")
H <- optimHess(f$opt$par, f$obj$fn, f$obj$gr)
D <- sqrt(abs(diag(H)))
ev <- eigen(H / outer(D, D), symmetric = TRUE)$values
cat("outer Hessian at stop: unit-diagonal eigenvalues", format(ev), "\n")
cat("hess_needs_reparam:", ns$hess_needs_reparam(H), "\n")
M0 <- ns$reparam_basis(f$frame, names(f$opt$par),
                       est = f$obj$env$parList(f$opt$par))
cat("reparam_basis NULL:", is.null(M0), "\n")
if (!is.null(M0)) print(M0)
X <- f$frame$linpreds[[1]]$X
cat("design_whitener NULL:", is.null(ns$design_whitener(X, 5e6)), "\n")
cat("predictor_design_kappa:", ns$predictor_design_kappa(f$frame), "\n")
cat("opt$polish:", format(unlist(f$opt$polish)), "\n")
# newton_polish from the stop, with a trace of its decisions
trace("hess_needs_reparam", where = ns, print = FALSE,
      exit = quote(cat("  hess_needs_reparam ->", returnValue(), "\n")))
trace("reparam_hessian_refined", where = ns, print = FALSE,
      exit = quote(cat("  reparam_hessian_refined NULL:",
                       is.null(returnValue()), "\n")))
op2 <- ns$newton_polish(f$obj, f$opt, list(lower = -Inf, upper = Inf),
                        f$par_units, f$control, ill = TRUE, frame = f$frame)
cat("polish again:", format(op2$par, digits = 10), "obj",
    format(op2$objective, digits = 12), "\n")
