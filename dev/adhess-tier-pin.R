# Lane adhess: how far the covariance of test-se-tier.R's "healthy fit
# with random effects" (sleepstudy, Reaction / 32, (Days | Subject))
# moves from 0.71.0's (optimHess() differences, absolute 1e-3) to the
# design-scaled steps (Days has a root mean square of 5.34), and both
# against a Richardson Jacobian of the AD gradient (numDeriv).
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
ss <- lme4::sleepstudy
ss$Reaction <- ss$Reaction / 32
fit <- frm(Reaction ~ Days + (Days | Subject), data = ss)
obj <- fit$obj
p <- ns$sdr_outer_point(obj)
st <- ns$obj_state_save(obj)
old <- solve(stats::optimHess(p, obj$fn, obj$gr))
ns$obj_state_restore(obj, st)
new <- unname(vcov(fit, full = TRUE))
H0 <- ns$fd_hessian(p, obj$gr)$H
D <- sqrt(abs(diag(H0)))
g0 <- function(dz) as.numeric(obj$gr(p + dz / D)) / D
J <- numDeriv::jacobian(g0, rep(0, length(p)), method = "Richardson",
                        method.args = list(eps = 0.01, d = 0.01,
                                           zero.tol = 1, r = 6, v = 2))
ref <- solve((J + t(J)) / 2 * outer(D, D))
ns$obj_state_restore(obj, st)
cat("move (old vs new):", format(cov_err(unname(old), new), digits = 3), "\n")
cat("old vs reference: ", format(cov_err(unname(old), ref), digits = 3), "\n")
cat("new vs reference: ", format(cov_err(new, ref), digits = 3), "\n")
