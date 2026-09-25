# As phase3b-cleft-probe-re.R, on the unfitted objective (dry_run).
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
obj <- if (!is.null(fit$obj)) fit$obj else fit; print(names(fit))
p0 <- obj$par
cat("par names:", names(p0), "\n")
cat("start:", format(p0, digits = 4), "\n")
cat("fn:", obj$fn(p0), " gr:", format(obj$gr(p0), digits = 4), "\n")
for (k in seq_along(p0)) for (h in c(-3, -1, -0.3, 0.3, 1, 3, 8)) {
  p <- p0; p[k] <- p[k] + h
  f <- obj$fn(p); g <- obj$gr(p)
  if (!is.finite(f) || any(!is.finite(g))) {
    cat("par", k, names(p0)[k], "step", h, "fn", f, "gr", format(g, digits = 3), "\n")
  }
}
env <- obj$env
full <- env$last.par
cat("full length", length(full), "random idx", length(env$random), "\n")
rnd <- env$random
for (sc in c(0, 0.5, 1, 2, 4, 8)) {
  for (k in c(0, 1)) {
    pp <- full; pp[rnd] <- if (k == 0) sc else -sc
    f <- env$f(pp, order = 0)
    g <- env$f(pp, order = 1)
    cat("u =", if (k == 0) sc else -sc, ": f", f, " grad finite", all(is.finite(g)), "\n")
  }
}
# which block
pp <- full; pp[rnd] <- 0
h <- tryCatch(env$spHess(pp, random = TRUE), error = function(e) conditionMessage(e))
cat("inner Hessian at u = 0 finite:", if (is.character(h)) h else all(is.finite(h@x)), "\n")
set.seed(1)
nbad <- 0
for (i in 1:40) {
  pp <- full; pp[rnd] <- rnorm(length(rnd), 0, c(1, 3, 10)[1 + i %% 3])
  h <- tryCatch(env$spHess(pp, random = TRUE), error = function(e) NULL)
  fin <- !is.null(h) && all(is.finite(h@x))
  g <- env$f(pp, order = 1)
  if (!fin || !all(is.finite(g))) {
    nbad <- nbad + 1
    if (nbad <= 3) cat("bad at draw", i, "u range", format(range(pp[rnd]), digits = 3),
                       "Hessian finite", fin, "grad finite", all(is.finite(g)), "\n")
  }
}
cat("bad draws", nbad, "of 40\n")
env$inner.control$trace <- 1
invisible(tryCatch(obj$fn(obj$par), error = function(e) cat("fn error", conditionMessage(e), "\n")))
cat("fn again:", obj$fn(obj$par), "\n")
cat("gr again:", obj$gr(obj$par), "\n")
cat("env$last.par.best random:", format(range(env$last.par.best[rnd]), digits = 3), "\n")
