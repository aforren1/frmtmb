# Reviewer probe: does wiener()'s identity-link support guard,
# 0 / (x + |x|), survive on the AD tape, where a fit evaluates it?
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("frmtmb from", find.package("frmtmb"), "; eam from",
    find.package("frmtmb.eam"), "\n")
cat("formals(wiener):", names(formals(wiener)), "\n")
f <- function(p) 0 / (p + abs(p))
tp <- RTMB::MakeTape(f, 1)
cat("tape of 0 / (p + |p|) at p = -1:", tp(-1), "; plain:", f(-1), "\n")

set.seed(20261016)
d <- frmtmb.eam:::ddm_simulate(60, mu = 0.8, bs = 1.3, ndt = 0.3)
fid <- frm(bf(rt | dec(upper) ~ 1),
           family = wiener(link_bs = "identity", link_ndt = "identity",
                           link_bias = "identity"), data = d)
cat("identity fit: logLik", as.numeric(logLik(fid)), "conv",
    fid$opt$convergence, "\n")
par <- fid$opt$par
cat("par names:", names(par), "\n")
cat("objective at the optimum:", fid$obj$fn(par), "\n")
for (nm in c("bias", "bs", "ndt")) {
  for (v in c(1.2, -0.5)) {
    p2 <- par
    i <- grep(paste0("^", nm), names(p2))
    if (!length(i)) i <- grep(nm, names(p2))
    if (!length(i)) next
    p2[i[1]] <- v
    cat(sprintf("objective (tape) with %s = %g: %s\n", nm, v,
                format(fid$obj$fn(p2))))
  }
}
# the plain-double density, as the lane's test checks it
y <- d$rt[1:3]
av <- list(dec = d$upper[1:3])
cat("plain lpdf at bias 1.2:",
    family(fid)$lpdf(y, list(mu = 0.8, bs = 1.3, ndt = 0.2, bias = 1.2), av),
    "\n")

cat("\n--- objective on the tape, betad moved outside the support\n")
print(fid$opt$par)
for (j in 1:3) {
  for (v in c(1.2, -0.5)) {
    p2 <- par
    p2[which(names(par) == "betad")[j]] <- v
    cat(sprintf("betad[%d] = %g: objective %s\n", j, v,
                format(fid$obj$fn(p2))))
  }
}
cat("\n--- data near the bias edge, identity link\n")
for (s in 1:5) {
  set.seed(20261100 + s)
  d2 <- frmtmb.eam:::ddm_simulate(80, mu = 0.5, bs = 1.0, ndt = 0.25,
                                  bias = 0.97)
  r <- tryCatch(withCallingHandlers(
    frm(bf(rt | dec(upper) ~ 1), family = wiener(link_bias = "identity"),
        data = d2),
    warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) e)
  if (inherits(r, "error")) { cat("seed", s, "ERROR", conditionMessage(r), "\n"); next }
  b <- r$opt$par[which(names(r$opt$par) == "betad")]
  cat(sprintf("seed %d: conv %d logLik %.4f betad %s\n", s,
              r$opt$convergence, as.numeric(logLik(r)),
              paste(signif(b, 5), collapse = " ")))
}
