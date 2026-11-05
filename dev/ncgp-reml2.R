# Why is the REML gradient large at low noise? Tape gradient against a central
# difference of the tape, by inner tolerance. Seed 20261102.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
s <- 1e-2
set.seed(20261102)
n0 <- 60L
x <- round(runif(n0, 0, 6), 3)
e0 <- rnorm(n0)
d <- data.frame(x = x, y = sin(2 * x) + s * e0)[1:40, ]
f <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian(), REML = TRUE))
fr <- f$frame
print(f$obj$env$inner.control)
p <- f$opt$par
tpl <- f$estimates
mk <- function(ctl) {
  frmtmb:::gp_tape_scope(fr, RTMB::MakeADFun(frmtmb:::build_objective(fr), tpl,
    random = c("b", "beta"), silent = TRUE, inner.control = ctl))
}
ng <- function(o, p, h = 1e-4) vapply(seq_along(p), function(i) {
  e <- replace(numeric(length(p)), i, h)
  (o$fn(p + e) - o$fn(p - e)) / (2 * h)
}, 0)
for (ctl in list(list(), list(grad.tol = 1e-12, step.tol = 1e-12), list(grad.tol = 1e-14, step.tol = 1e-14, tol10 = 1e-8))) {
  o <- mk(ctl)
  cat("ctl", unlist(ctl), ": fn", format(o$fn(p), digits = 12), " gr", format(o$gr(p), digits = 4), " fd", format(ng(o, p), digits = 4), "\n")
  op <- nlminb(p, o$fn, o$gr, control = list(iter.max = 200))
  cat("   nlminb:", op$message, format(op$objective, digits = 12), format(op$par, digits = 6), "\n")
}
# the ML tape at the same data: gradient against fd
f1 <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
o <- f1$obj; p1 <- f1$opt$par
cat("ML: gr", format(o$gr(p1), digits = 4), " fd", format(ng(o, p1), digits = 4), "\n")
inner <- o$env$last.par[names(o$env$last.par) == "b"]
cat("ML inner gradient max", max(abs(o$env$f(o$env$last.par, order = 1)[names(o$env$last.par) == "b"])), "\n")
o2 <- f$obj; o2$fn(p)
lp <- o2$env$last.par
gi <- o2$env$f(lp, order = 1)
cat("REML inner gradient: b", max(abs(gi[names(lp) == "b"])), "beta", max(abs(gi[names(lp) == "beta"])), "\n")
h <- o2$env$spHess(lp, random = TRUE)
ev <- eigen(as.matrix(h), only.values = TRUE)$values
cat("REML inner H cond", max(ev) / min(ev), "\n")
cat("--- arms\n")
nll <- frmtmb:::build_objective(fr)
arm <- function(comp, atomic_all) {
  old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  tc <- RTMB::TapeConfig()
  on.exit({TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB"); RTMB::TapeConfig(tc)})
  TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
  if (atomic_all) RTMB::TapeConfig(matmul = "atomic")
  o <- RTMB::MakeADFun(nll, tpl, random = c("b", "beta"), silent = TRUE)
  force(o$env$spHess)
  o
}
for (comp in 0:1) for (aa in c(FALSE, TRUE)) {
  o <- arm(comp, aa)
  cat("compress", comp, "atomic_all", aa, ": gr", format(o$gr(p), digits = 5), " fd", format(ng(o, p), digits = 5), "\n")
}
