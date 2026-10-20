# Lane adhess, finding 8: where the ML fit on Days + 1e5 stops, and what
# newton_polish() does there (verbose fit), against lme4's maximum.
#   ADHESS_ARM=lane|base Rscript dev/adhess-shiftdbg.R <shift>
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
sh <- as.numeric(commandArgs(TRUE)[1])
data("sleepstudy", package = "lme4")
set.seed(1)
d <- sleepstudy[sort(sample(nrow(sleepstudy), 120)), ]
d$t <- d$Days + sh
f <- suppressWarnings(frm(Reaction ~ t + (1 | Subject), data = d,
                          verbose = TRUE))
print(f$opt[c("par", "objective", "convergence", "message", "polish")])
l <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = d, REML = FALSE)
cat("lme4 -logLik", -as.numeric(logLik(l)), "\n")
p <- f$opt$par
g <- f$obj$gr(p)
cat("gradient at the stop:", format(g, digits = 4), "\n")
cat("design kappa", ns$predictor_design_kappa(f$frame), "\n")
cat("par_units:", format(f$par_units %||% 1, digits = 4), "\n")
M <- ns$reparam_basis(f$frame, names(p))
print(M)
h <- ns$reparam_hessian_refined(f$obj, p, M)
st <- -solve(h$H, drop(crossprod(h$M, as.numeric(g))))
cat("Newton gain in the design's coordinates:",
    -sum(drop(crossprod(h$M, as.numeric(g))) * st) / 2, "\n")
q <- p + drop(h$M %*% st)
cat("objective after that step:", f$obj$fn(q), "\n")
# newton_polish() by hand from the stop, with its Hessian decisions shown
trace("hess_needs_reparam", exit = quote(cat("hess_needs_reparam ->",
                                             returnValue(), "\n")),
      where = ns, print = FALSE)
trace("reparam_hessian_refined", exit = quote({
  rv <- returnValue()
  cat("reparam H:", if (is.null(rv)) "NULL" else
        format(eigen(rv$H, only.values = TRUE)$values, digits = 3), "\n")
}), where = ns, print = FALSE)
op <- ns$newton_polish(f$obj, f$opt, f$cache$bounds, f$par_units,
                       f$control, ill = TRUE, frame = f$frame)
cat("polish:", format(op$polish), "objective", op$objective, "\n")
