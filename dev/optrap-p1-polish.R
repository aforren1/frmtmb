# Lane optrap, punch round 1 (m2, m10): a candidate, not in the package.
# After the fit, Newton steps on TMB's exact Hessian (fixed effects
# only, no random effects here), each kept only if the objective falls,
# up to 5. Raw polynomials of dev/optrap-p1-polystudy.R, degrees 5 and
# 6, seeds 1 to 20: the gap to lm() before and after, and the largest
# |beta - lm| / SE after.
#   Rscript dev/optrap-p1-polish.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
polish <- function(f) {
  p <- f$opt$par
  v <- f$obj$fn(p)
  u <- f$par_units %||% rep(1, length(p))
  for (it in 1:5) {
    g <- drop(f$obj$gr(p)) * u
    H <- f$obj$he(p) * outer(u, u)
    st <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(st) || !all(is.finite(st))) break
    q <- p + st * u
    vq <- f$obj$fn(q)
    if (!is.finite(vq) || vq >= v) break
    p <- q
    v <- vq
  }
  c(obj = v, steps = it)
}
for (deg in 5:6) {
  r <- NULL
  for (s in 1:20) {
    set.seed(1000 * deg + s)
    d <- data.frame(x = runif(200, 1, 2))
    d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
    for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
    fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
    f <- suppressWarnings(frm(fo, data = d))
    l <- as.numeric(logLik(lm(fo, data = d)))
    pl <- polish(f)
    r <- rbind(r, c(before = -f$opt$objective - l, after = -pl[["obj"]] - l))
  }
  cat(sprintf("deg %d | short before %2d (min %8.3g) | short after %2d (min %8.3g)\n",
              deg, sum(r[, 1] < -1e-4), min(r[, 1]), sum(r[, 2] < -1e-4),
              min(r[, 2])))
}
