# Item 3.4: validate the Wiener response-time distribution function.
#
#  1. Against the 400-bit reference (dev/phase3b-cdf-reference-v2.R), on
#     the grid the density is pinned on and a tail extension: relative
#     error of F and of S, each on the log scale the family returns.
#  2. The blend center sweep that chose ddm_rtcdf_u0.
#  3. Against RWiener::pwiener(resp = "both") and WienR::pWDM(), the two
#     established implementations, in absolute terms (their own
#     accuracy is absolute), and the reference's verdict on THEM.
#  4. F integrates the density: F(t2) - F(t1) against integrate() of
#     the family's own density summed over both boundaries.
#
# Output: dev/phase3b-log/cdf-validate.txt
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (nzchar(Sys.getenv("P3B_INSTALLED"))) {
  ns <- asNamespace("frmtmb.eam")
  for (nm in ls(ns, all.names = TRUE)) assign(nm, get(nm, ns))
} else {
  suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE,
                                     export_all = TRUE))
}
out <- character(0)
say <- function(...) {
  s <- paste0(...)
  out <<- c(out, s)
  cat(s, "\n")
}

ref <- utils::read.csv("dev/phase3b-log/cdf-reference-400.csv",
                       colClasses = c(F_ref = "character",
                                      S_ref = "character"))
Fr <- as.numeric(ref$F_ref)
Sr <- as.numeric(ref$S_ref)
u <- ref$t / ref$a^2
say("rows: ", nrow(ref), " (", sum(ref$grid == "density"), " density grid, ",
    sum(ref$grid == "tail"), " tail); u from ", signif(min(u), 3), " to ",
    signif(max(u), 3))
say("F from ", signif(min(Fr), 3), "; S from ", signif(min(Sr), 3))

blend <- function(u0, us = 0.12) {
  tt <- ref$t
  uu <- tt / ref$a^2
  lam <- 0.5 * (1 + tanh((log(uu) - log(u0)) / us))
  fs <- ddm_cdf_small(tt, ref$v, ref$a, ref$w)
  sl <- ddm_surv_large(tt, ref$v, ref$a, ref$w)
  lF <- (1 - lam) * log(ddm_floor(fs, ddm_share_floor)) +
    lam * log(ddm_floor(1 - sl, ddm_share_floor))
  lS <- (1 - lam) * log(ddm_floor(1 - fs, ddm_share_floor)) +
    lam * log(ddm_floor(sl, ddm_share_floor))
  list(F = exp(lF), S = exp(lS))
}
rel <- function(x, r) abs(x - r) / r

say("")
say("== blend center sweep, max relative error over all rows ==")
say("(one center for both columns here; the shipped function takes F at")
say(" ddm_rtcdf_u0F and S at ddm_rtcdf_u0S, each from its own column)")
say(sprintf("%6s %12s %12s %8s %8s", "u0", "F max rel", "S max rel",
            "F>1e-10", "S>1e-10"))
for (u0 in c(0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.35)) {
  b <- blend(u0)
  eF <- rel(b$F, Fr)
  eS <- rel(b$S, Sr)
  say(sprintf("%6.3f %12.3e %12.3e %8d %8d", u0, max(eF), max(eS),
              sum(eF > 1e-10), sum(eS > 1e-10)))
}

say("")
say("== the shipped function, ddm_rt_lcdf2() ==")
r <- ddm_rt_lcdf2(ref$t, ref$v, ref$a, ref$w)
eF <- rel(exp(r$lF), Fr)
eS <- rel(exp(r$lS), Sr)
for (g in c("density", "tail")) {
  k <- ref$grid == g
  say(g, " grid: F max rel ", signif(max(eF[k]), 3), ", S max rel ",
      signif(max(eS[k]), 3), ", F max abs ",
      signif(max(abs(exp(r$lF[k]) - Fr[k])), 3))
}
w <- which.max(eF)
say("worst F row: t=", ref$t[w], " v=", ref$v[w], " a=", ref$a[w], " w=",
    ref$w[w], " F=", signif(Fr[w], 4))
w <- which.max(eS)
say("worst S row: t=", ref$t[w], " v=", ref$v[w], " a=", ref$a[w], " w=",
    ref$w[w], " S=", signif(Sr[w], 4))

say("")
say("== established implementations, against the same 400-bit truth ==")
rw <- mapply(function(t, v, a, w) {
  RWiener::pwiener(t + 1e-9, a, 1e-9, w, v, resp = "both")
}, ref$t, ref$v, ref$a, ref$w)
wr <- mapply(function(t, v, a, w) {
  WienR::pWDM(t, "lower", a = a, v = v, w = w, precision = 1e-12)$value +
    WienR::pWDM(t, "upper", a = a, v = v, w = w, precision = 1e-12)$value
}, ref$t, ref$v, ref$a, ref$w)
kd <- ref$grid == "density"
say("frmtmb.eam vs RWiener::pwiener, density grid, max abs: ",
    signif(max(abs(exp(r$lF[kd]) - rw[kd])), 3))
say("frmtmb.eam vs RWiener::pwiener, all rows, max abs: ",
    signif(max(abs(exp(r$lF) - rw)), 3))
say("RWiener vs the 400-bit truth, max abs: ", signif(max(abs(rw - Fr)), 3),
    "; max rel: ", signif(max(rel(rw, Fr)), 3))
say("frmtmb.eam vs WienR::pWDM(precision = 1e-12), all rows, max abs: ",
    signif(max(abs(exp(r$lF) - wr)), 3))
say("WienR vs the 400-bit truth, max abs: ", signif(max(abs(wr - Fr)), 3),
    "; max rel: ", signif(max(rel(wr, Fr)), 3))
say("frmtmb.eam vs the 400-bit truth, max abs: ",
    signif(max(abs(exp(r$lF) - Fr)), 3))

say("")
say("== F integrates the density ==")
worst <- 0
set.seed(34)
for (i in 1:40) {
  v <- runif(1, -3, 3); a <- runif(1, 0.5, 3); w <- runif(1, 0.2, 0.8)
  t1 <- runif(1, 0.01, 0.5); t2 <- t1 + runif(1, 0.05, 3)
  dens <- function(t) {
    exp(ddm_lpdf_both(t, v, a, w, 0)) + exp(ddm_lpdf_both(t, v, a, w, 1))
  }
  I <- stats::integrate(dens, t1, t2, rel.tol = 1e-12,
                        subdivisions = 2000L)$value
  Fd <- diff(exp(ddm_rt_lcdf2(c(t1, t2), v, a, w)$lF))
  worst <- max(worst, abs(Fd - I) / I)
}
say("40 random (v, a, w, t1, t2): max |F(t2) - F(t1) - integral| / ",
    "integral = ", signif(worst, 3))
writeLines(out, "dev/phase3b-log/cdf-validate.txt")
