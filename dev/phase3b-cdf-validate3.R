# Punch round 1, M1 and item (b): validate the rewritten response-time
# distribution against the 700-bit reference of
# dev/phase3b-cdf-reference-v3.R, over v in -20..20, a in 0.3..6,
# w in 0.001..0.999, u in 1e-3..10 (3150 rows).
#
#  1. Blend-center sweeps for S, F and the defective F_b, error by |v| a.
#  2. The shipped functions, worst error by |v| a band, in log space
#     (absolute error of the log, which is the relative error of the
#     quantity where that is small).
#  3. F_lower + F_upper = F, the per-boundary functions against the
#     marginal one, at double precision.
#  4. RWiener::pwiener() and WienR::pWDM() per boundary, against the
#     same reference.
# Output: dev/phase3b-log/cdf-validate3.txt
.libPaths(c(Sys.getenv("P3B_LIB", "C:/Users/adf44/source/r/phase3b-lib"),
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.eam")
out <- character(0)
say <- function(...) {
  s <- paste0(...)
  out <<- c(out, s)
  cat(s, "\n")
}
ref <- utils::read.csv("dev/phase3b-log/cdf-reference-v3.csv",
                       colClasses = c(lS = "character", lF = "character",
                                      lFl = "character", lFu = "character"))
for (nm in c("lS", "lF", "lFl", "lFu")) ref[[nm]] <- as.numeric(ref[[nm]])
say("package: ", find.package("frmtmb.eam"))
say("rows: ", nrow(ref))
va <- abs(ref$v) * ref$a
band <- cut(va, c(-1, 1, 5, 12, 24, 48, 72, 121),
            labels = c("<=1", "<=5", "<=12", "<=24", "<=48", "<=72", "<=120"))
t <- as.numeric(ref$t); v <- as.numeric(ref$v); a <- as.numeric(ref$a)
w <- as.numeric(ref$w)
lu <- log(ref$u)

# log-scale error; for a quantity whose log is near zero (F near 1) the
# absolute log error IS the relative error, and elsewhere it is the
# relative error too, to first order
lerr <- function(x, r) abs(x - r)

# the routes, once
Ss <- ns$ddm_lsurv_small(t, v, a, w)
Sl <- ns$ddm_lsurv_large(t, v, a, w)
Fls <- ns$ddm_llower_small(t, v, a, w)
Fus <- ns$ddm_llower_small(t, -v, a, 1 - w)
Fll <- ns$ddm_llower_large(t, v, a, w)
Ful <- ns$ddm_llower_large(t, -v, a, 1 - w)
Fs <- ns$ddm_lse2(Fls, Fus)
FlS <- log(pmax(-expm1(pmin(Sl, 0)), 1e-300))
lam <- function(u0) 0.5 * (1 + tanh((lu - log(u0)) / ns$ddm_rtcdf_us))

say("")
say("== blend-center sweep: worst log error over all 3150 rows ==")
say(sprintf("%6s %10s %10s %10s %10s", "u0", "S", "F", "F_lower", "F_upper"))
for (u0 in c(0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 1)) {
  L <- lam(u0)
  say(sprintf("%6.2f %10.3g %10.3g %10.3g %10.3g", u0,
              max(lerr((1 - L) * Ss + L * Sl, ref$lS)),
              max(lerr((1 - L) * Fs + L * FlS, ref$lF)),
              max(lerr((1 - L) * Fls + L * Fll, ref$lFl)),
              max(lerr((1 - L) * Fus + L * Ful, ref$lFu))))
}

say("")
say(sprintf("== shipped: u0F = %g, u0S = %g, u0B = %g ==", ns$ddm_rtcdf_u0F,
            ns$ddm_rtcdf_u0S, ns$ddm_rtcdf_u0B))
r <- ns$ddm_rt_lcdf2(t, v, a, w)
bl <- ns$ddm_rt_lcdf_b(t, v, a, w, 0)
bu <- ns$ddm_rt_lcdf_b(t, v, a, w, 1)
e <- data.frame(band, S = lerr(r$lS, ref$lS), F = lerr(r$lF, ref$lF),
                Fl = lerr(bl, ref$lFl), Fu = lerr(bu, ref$lFu))
tab <- aggregate(cbind(S, F, Fl, Fu) ~ band, data = e, FUN = max)
tab$rows <- as.integer(table(band)[as.character(tab$band)])
say(paste(utils::capture.output(print(tab, digits = 3, row.names = FALSE)),
          collapse = "\n"))
for (nm in c("S", "F", "Fl", "Fu")) {
  i <- which.max(e[[nm]])
  say(sprintf("worst %s: v=%g a=%g w=%g u=%g, reference log %.6g, got %.6g",
              nm, v[i], a[i], w[i], ref$u[i],
              ref[[c(S = "lS", F = "lF", Fl = "lFl", Fu = "lFu")[[nm]]]][i],
              c(S = r$lS[i], F = r$lF[i], Fl = bl[i], Fu = bu[i])[[nm]]))
}
say(sprintf("smallest reference log S reproduced: %.1f (got %.1f)",
            min(ref$lS), r$lS[which.min(ref$lS)]))

say("")
say("== the per-boundary functions sum to the marginal one ==")
sm <- ns$ddm_lse2(bl, bu)
say("max |log(F_lower + F_upper) - log F|: ", format(max(abs(sm - r$lF)),
                                                    digits = 3))

say("")
say("== established implementations per boundary, against the reference ==")
k <- ref$u >= 5e-3
rwl <- mapply(function(t, v, a, w) {
  RWiener::pwiener(t + 1e-9, a, 1e-9, w, v, resp = "lower")
}, t[k], v[k], a[k], w[k])
wrl <- mapply(function(t, v, a, w) {
  WienR::pWDM(t, "lower", a = a, v = v, w = w, precision = 1e-12)$value
}, t[k], v[k], a[k], w[k])
say(sprintf("rows with u >= 5e-3: %d", sum(k)))
say("RWiener lower, max abs error on F_lower: ",
    format(max(abs(rwl - exp(ref$lFl[k]))), digits = 3),
    "; frmtmb.eam: ", format(max(abs(exp(bl[k]) - exp(ref$lFl[k]))), digits = 3))
say("WienR lower (precision = 1e-12), max abs error on F_lower: ",
    format(max(abs(wrl - exp(ref$lFl[k]))), digits = 3),
    "; frmtmb.eam max log error on the same rows: ",
    format(max(abs(bl[k] - ref$lFl[k])), digits = 3))
writeLines(out, "dev/phase3b-log/cdf-validate3.txt")

say("")
say("== tape gradients against central differences ==")
set.seed(35)
worst <- 0
for (i in 1:60) {
  p <- c(runif(1, -8, 8), log(runif(1, 0.4, 4)), runif(1, 0.02, 0.2))
  q <- c(0.3, 0.6, 1.5); up <- c(0, 1, 1)
  for (nm in c("lS", "lF", "lFb")) {
    f <- function(p) {
      if (nm == "lFb") return(sum(ns$ddm_rt_lcdf_b(q - p[3], p[1], exp(p[2]),
                                                    0.4, up)))
      sum(ns$ddm_rt_lcdf2(q - p[3], p[1], exp(p[2]), 0.4)[[nm]])
    }
    tp <- RTMB::MakeTape(f, p)
    g <- tp$jacobian(p)
    fd <- vapply(1:3, function(k) {
      h <- 1e-6 * max(1, abs(p[k])); e <- replace(numeric(3), k, h)
      (f(p + e) - f(p - e)) / (2 * h)
    }, 0)
    worst <- max(worst, max(abs(g - fd) / pmax(abs(fd), 1)))
  }
}
say("60 random points, three functions: max |tape - central difference| / max(|fd|, 1) = ",
    format(worst, digits = 3))
writeLines(out, "dev/phase3b-log/cdf-validate3.txt")
