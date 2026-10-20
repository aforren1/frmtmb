# Review of lane illse2, item 4: an interior maximum just inside a bound
# must not be pinned to it. The bound sits c standard errors beyond
# lm()'s estimate of the top coefficient of a raw polynomial, so the
# unconstrained maximum (lm()) is the constrained one.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-bound.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
capm <- function(expr) {
  w <- character(); ms <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "i2_err")),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") },
    message = function(x) { ms <<- c(ms, conditionMessage(x))
      invokeRestart("muffleMessage") })
  list(v = v, w = w, m = ms)
}
for (deg in c(2L, 6L, 7L)) for (seed in 1:3) for (cc in c(1e-3, 1e-2, 0.1)) {
  for (side in c("ub", "lb")) {
    set.seed(seed)
    d <- data.frame(x = runif(200, 1, 2))
    d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
    xn <- paste0("x", 1:deg)
    for (k in 1:deg) d[[xn[k]]] <- d$x^k
    top <- xn[deg]
    fo <- as.formula(paste("y ~", paste0(xn, collapse = " + ")))
    l0 <- lm(fo, data = d)
    b <- coef(l0)[[top]]
    s <- sqrt(diag(vcov(l0)))[[top]]
    n <- nrow(d)
    se_ml <- s * sqrt((n - deg - 1) / n)
    pr <- if (side == "ub") {
      set_prior("", class = "b", coef = top, ub = b + cc * s)
    } else set_prior("", class = "b", coef = top, lb = b - cc * s)
    r <- capm(frm(fo, data = d, prior = pr))
    if (is_err(r$v)) { cat(deg, seed, cc, side, "ERROR", r$v, "\n"); next }
    f <- r$v
    fe <- suppressWarnings(fixef(f))
    lost <- ns$sdr_of(f)$se_lost
    cat(sprintf(paste0("deg %d seed %d c %-5g %s | code %d | (est - lm)/SE %+.2e | ",
                       "logLik - lm %+.2e | top SE / lm %.6f | lost %s | ",
                       "warn %s | msg %s\n"),
                deg, seed, cc, side, f$opt$convergence,
                (fe[top, "Estimate"] - b) / s,
                as.numeric(logLik(f)) - as.numeric(logLik(l0)),
                fe[top, "Est.Error"] / se_ml,
                if (length(lost)) paste(names(lost), lost, collapse = ",")
                else "none", wtag(r$w),
                if (length(r$m)) substr(paste(r$m, collapse = " | "), 1, 90)
                else "-"))
  }
}
