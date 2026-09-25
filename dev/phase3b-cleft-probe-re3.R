# The FIRST objective call on the unfitted object, traced.
# steps, then evaluate the objective and gradient at the start and along
# each coordinate.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
set.seed(1)
NS <- 8L; NT <- 400L
u <- rnorm(NS, 0, 0.35); b <- rnorm(NS, 0, 0.2)
d <- do.call(rbind, lapply(seq_len(NS), function(s) {
  cond <- rep(0:1, length.out = NT)
  x <- ddm_simulate(NT, mu = 0.4 + 0.9 * cond + u[s],
                    bs = 1.4 * exp(b[s]), ndt = 0.25, bias = 0.5)
  x$cond <- cond; x$s <- factor(s); x
}))
d$code <- 0L; d$y2 <- d$rt
fast <- d$rt < 0.45
d$code[fast] <- -1L; d$rt[fast] <- 0.45
pick <- which(!fast)[seq(5, sum(!fast), by = 5)]
lo <- floor(d$rt[pick] * 10) / 10
d$code[pick] <- 2L; d$y2[pick] <- lo + 0.1; d$rt[pick] <- pmax(lo, 0.45)
fit <- suppressWarnings(frm(bf(rt | dec(upper) + cens(code, y2) ~ cond + (1 | s),
                               bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
                            family = wiener(), data = d, dry_run = "objective",
                            control = frmtmb_control(optCtrl = list(
                              iter.max = 0, eval.max = 1), restarts = 0)))
obj <- fit$obj
env <- obj$env
cat("initial random:", format(range(env$last.par[env$random]), digits = 3),
    " par:", format(env$last.par[-env$random], digits = 3), "\n")
cat("f at initial full:", env$f(env$last.par, order = 0), "\n")
env$inner.control$trace <- 1
cat("fn first call:", obj$fn(obj$par), "\n")
full <- env$last.par; full[env$random] <- 0
g <- env$f(full, order = 1)[env$random]
H <- as.matrix(env$spHess(full, random = TRUE))
step <- -solve(H, g)
cat("Newton step range:", format(range(step), digits = 3), "\n")
for (s in c(1, 0.5, 0.1, 0.01, 1e-3)) {
  pp <- full; pp[env$random] <- s * step
  cat("step x", s, ": f", env$f(pp, order = 0), "\n")
}
for (k in seq_along(env$random)) for (val in c(-50, -20, 20, 50)) {
  pp <- full; pp[env$random[k]] <- val
  f <- env$f(pp, order = 0)
  if (!is.finite(f)) cat("u[", k, "] =", val, "f", f, "\n")
}
for (s in c(1, 0.5, 0.1)) {
  pp <- full; pp[env$random] <- s * step
  g <- env$f(pp, order = 1)
  h <- tryCatch(as.matrix(env$spHess(pp, random = TRUE)), error = function(e) NULL)
  cat("step x", s, ": grad finite", all(is.finite(g)), " Hessian finite",
      !is.null(h) && all(is.finite(h)), "\n")
}
g0 <- env$f(full, order = 1)
cat("grad at 0 finite", all(is.finite(g0)), " H at 0 finite", all(is.finite(H)), "\n")
print(str(env$inner.control))
ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
cat("eigen of inner H at 0: min", min(ev), "max", max(ev), "\n")
# which diagonal entries
cat("diag H:", format(diag(H), digits = 3), "\n")
